"""The verifier driver, and the cost model checked against the verifiers' own numbers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from quantum_compilation.cost import DEFAULT_CALIBRATION, MODEL_AGREEMENT_REL, TASKS, score
from quantum_compilation.submission import audit
from quantum_compilation.verify import (
    VERIFIER_SCRIPTS,
    Verification,
    candidate_hashes,
    measured_calibration,
    model_agrees,
    parse_verifier_output,
    passing_costs,
    recompute_cost,
    run_verifier,
    verifier_script,
    verify_submission,
    write_rows,
)

REFERENCE_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def reference() -> dict[str, Verification]:
    """The three official verdicts on `sample_submission/`, replayed through the parser."""
    return {
        task: parse_verifier_output(
            task, (REFERENCE_DIR / f"reference_{task}.json").read_text(), 0.0
        )
        for task in TASKS
    }


def fake_verifiers(root: Path, task: str, body: str) -> Path:
    """A verifiers tree whose `task` script runs `body`, for testing the driver itself."""
    script = root / VERIFIER_SCRIPTS[task]
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(body)
    return root


class TestModelAgainstTheOfficialVerifiers:
    """`cost` reproduces each verifier's own score from the fields that verifier reported.

    The load-bearing test here: it is what distinguishes a transcription of the Evaluation
    page from the metric the leaderboard actually runs. The page omits the rotation price,
    that the Toffoli count is a simulated mean, and that the dqi failure rate is a
    confidence bound -- each of which these captured numbers would expose.
    """

    @pytest.mark.parametrize("task", TASKS)
    def test_recomputes_the_verifier_score(self, task, reference):
        verdict = reference[task]
        assert verdict.passed
        assert recompute_cost(verdict) == pytest.approx(
            verdict.score, rel=MODEL_AGREEMENT_REL
        )
        assert model_agrees(verdict) is True

    def test_the_reference_submission_scores_exactly_one(self, reference):
        # True by definition -- which makes it a check on the calibration, not the metric.
        costs = passing_costs(list(reference.values()))
        assert costs == pytest.approx(DEFAULT_CALIBRATION.reference_costs)
        assert score(costs, DEFAULT_CALIBRATION) == pytest.approx(1.0)

    def test_the_builtin_calibration_names_the_circuits_it_measured(self, reference):
        assert DEFAULT_CALIBRATION.reference_sha256 == candidate_hashes(
            list(reference.values())
        )

    @pytest.mark.parametrize(
        ("task", "strictly_below"), [("invsqrt", False), ("dqi", True)]
    )
    def test_the_static_toffoli_count_only_bounds_the_scored_one(
        self, task, strictly_below, reference
    ):
        # Why `resources` cannot price a circuit: measurement-based uncomputation makes the
        # executed count data-dependent, so the static count is an upper bound.
        detail = reference[task].detail
        assert detail["toffoli_count"] <= detail["max_toffoli_count"]
        if strictly_below:
            assert detail["toffoli_count"] < detail["max_toffoli_count"]


class TestOutputParsing:
    """The verifiers' shared JSON contract, and what a driver must do when it is not met."""

    def test_a_failure_is_not_a_pass_and_carries_its_reason(self):
        verdict = parse_verifier_output(
            "dqi", json.dumps({"status": "FAIL", "error": "eps_plus >= 1/16"}), 1.0
        )
        assert not verdict.passed
        assert verdict.error == "eps_plus >= 1/16"

    def test_a_pass_without_a_score_is_not_a_pass(self):
        # Treating a scoreless PASS as scored would put None into a cost calculation.
        assert not parse_verifier_output("dqi", json.dumps({"status": "PASS"}), 1.0).passed

    @pytest.mark.parametrize("stdout", ["", "Killed", "[1, 2]", "{truncated"])
    def test_output_that_is_not_a_json_object_is_unparseable(self, stdout):
        verdict = parse_verifier_output("invsqrt", stdout, 1.0)
        assert verdict.status == "UNPARSEABLE"
        assert verdict.score is None

    def test_an_integer_score_is_accepted_as_a_cost(self):
        verdict = parse_verifier_output(
            "dqi", json.dumps({"status": "PASS", "score": 7}), 1.0
        )
        assert verdict.passed
        assert verdict.score == 7.0

    @pytest.mark.parametrize("key", ["circuit_sha256", "submission_sha256"])
    def test_the_candidate_hash_is_read_under_either_name(self, key):
        # invsqrt and dqi report one, dfthcsos the other.
        verdict = parse_verifier_output(
            "dqi", json.dumps({"status": "PASS", "score": 1.0, key: "abc"}), 1.0
        )
        assert verdict.candidate_sha256 == "abc"

    def test_a_model_that_cannot_be_evaluated_is_not_a_disagreement(self):
        # A PASS whose detail lacks the formula's inputs must read as "not evaluated", not
        # as "the model is wrong" and not as agreement.
        verdict = parse_verifier_output(
            "dqi", json.dumps({"status": "PASS", "score": 5.0}), 1.0
        )
        assert recompute_cost(verdict) is None
        assert model_agrees(verdict) is None


class TestDriver:
    """Locating each verifier, pointing it at the right target, and running it."""

    def test_each_task_has_a_script_path(self):
        assert set(VERIFIER_SCRIPTS) == set(TASKS)

    def test_a_missing_verifiers_tree_says_what_to_point_at(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="verifiers"):
            verifier_script(tmp_path, "invsqrt")

    def test_an_unknown_task_is_rejected(self, tmp_path):
        with pytest.raises(KeyError, match="unknown task"):
            verifier_script(tmp_path, "nope")

    def test_a_target_the_audit_did_not_find_is_reported_as_missing(self, tmp_path):
        verdicts = list(verify_submission(tmp_path, audit(tmp_path)))
        assert [v.status for v in verdicts] == ["MISSING"] * len(TASKS)

    @pytest.mark.parametrize(
        ("body", "status"),
        [
            ("import json; print(json.dumps({'status': 'PASS', 'score': 12.5}))", "PASS"),
            ("raise SystemExit('no kickmix')", "CRASH"),
        ],
    )
    def test_runs_the_verifier_and_reports_how_it_ended(self, tmp_path, body, status):
        root = fake_verifiers(tmp_path / "verifiers", "invsqrt", body)
        verdict = run_verifier(root, "invsqrt", tmp_path / "invsqrt.kmb")
        assert verdict.status == status

    def test_a_crash_keeps_the_reason_from_stderr(self, tmp_path):
        root = fake_verifiers(tmp_path / "verifiers", "dqi", "raise SystemExit('no kickmix')")
        assert "no kickmix" in run_verifier(root, "dqi", tmp_path / "dqi.kmb").error

    def test_a_run_that_overruns_its_timeout_is_reported_not_awaited(self, tmp_path):
        root = fake_verifiers(tmp_path / "verifiers", "dqi", "import time; time.sleep(30)")
        verdict = run_verifier(root, "dqi", tmp_path / "dqi.kmb", timeout=0.5)
        assert verdict.status == "TIMEOUT"
        assert verdict.score is None

    def test_the_salt_and_seed_reach_the_verifier_and_are_recorded(self, tmp_path):
        # The seed is the comparability mechanism: all three verifiers derive their own from
        # the candidate's hash unless one is passed, so a sweep must be able to pin it.
        root = fake_verifiers(
            tmp_path / "verifiers",
            "invsqrt",
            "import json, sys; print(json.dumps({'status': 'PASS', 'score': 1.0, "
            "'argv': sys.argv[1:]}))",
        )
        verdict = run_verifier(
            root, "invsqrt", tmp_path / "c.kmb", salt="dev-7", seed=1234
        )
        assert verdict.detail["argv"][-4:] == ["--salt", "dev-7", "--seed", "1234"]
        assert (verdict.salt, verdict.seed) == ("dev-7", 1234)

    def test_no_seed_flag_is_passed_when_none_is_asked_for(self, tmp_path):
        root = fake_verifiers(
            tmp_path / "verifiers",
            "invsqrt",
            "import json, sys; print(json.dumps({'status': 'PASS', 'score': 1.0, "
            "'argv': sys.argv[1:]}))",
        )
        verdict = run_verifier(root, "invsqrt", tmp_path / "c.kmb")
        assert "--seed" not in verdict.detail["argv"]


class TestResultRows:
    """Persisting verdicts, which is what makes two runs comparable after the fact."""

    def test_a_row_records_the_invocation_not_just_the_score(self, reference):
        # A cost is a simulated mean at a seed; a row without the seed cannot be compared.
        row = reference["dqi"].row()
        assert row["score"] == pytest.approx(reference["dqi"].score)
        assert row["model_cost"] == pytest.approx(row["score"], rel=MODEL_AGREEMENT_REL)
        assert set(row) >= {"task", "status", "salt", "seed", "candidate_sha256", "seconds"}

    def test_rows_are_rewritten_after_each_verdict_so_an_interrupt_keeps_what_ran(
        self, tmp_path, reference
    ):
        path = tmp_path / "rows.jsonl"
        landed = [reference["invsqrt"]]
        write_rows(path, landed)
        assert len(path.read_text().strip().splitlines()) == 1
        write_rows(path, landed + [reference["dqi"]])
        assert len(path.read_text().strip().splitlines()) == 2


class TestMeasuredCalibration:
    """Turning three passing reference verdicts into the calibration -- gate G0's output."""

    def test_passing_verdicts_become_a_measured_calibration(self, reference):
        calibration = measured_calibration(list(reference.values()), "test")
        assert calibration.measured is True
        assert calibration.reference_costs == pytest.approx(
            DEFAULT_CALIBRATION.reference_costs
        )

    def test_the_calibration_records_which_circuits_produced_it(self, reference):
        # Without this, calibrating against a compressed candidate would silently make that
        # candidate score 1.0 and every later ratio relative to it.
        calibration = measured_calibration(list(reference.values()), "test")
        assert calibration.reference_sha256 == DEFAULT_CALIBRATION.reference_sha256
        assert calibration.disagreeing_hashes({"dqi": "something-else"}) == {
            "dqi": (DEFAULT_CALIBRATION.reference_sha256["dqi"], "something-else")
        }

    def test_refuses_to_calibrate_from_an_incomplete_run(self, reference):
        failed = Verification("dqi", "FAIL", None, 1.0, detail={"error": "rejected"})
        partial = [reference["invsqrt"], failed, reference["dfthcsos"]]
        with pytest.raises(ValueError, match="must pass to calibrate"):
            measured_calibration(partial, "test")
