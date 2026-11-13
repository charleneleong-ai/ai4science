"""The verifier driver, and the cost model checked against the verifiers' own numbers."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from quantum_compilation.cost import DEFAULT_CALIBRATION, MODEL_AGREEMENT_REL, TASKS, score
from quantum_compilation.submission import SubmissionAudit
from quantum_compilation.verify import (
    VERIFIER_SCRIPTS,
    Verification,
    VerificationRun,
    load_verifier,
    parse_verifier_output,
    run_verifier,
    verifier_script,
)

REFERENCE_DIR = Path(__file__).parent / "data"

# A verifier body that passes and reports back the argv it was handed, for the two tests that
# check what the driver puts on the command line.
ECHO_ARGV = (
    "import json, sys; print(json.dumps({'status': 'PASS', 'score': 1.0, "
    "'argv': sys.argv[1:]}))"
)


@pytest.fixture(scope="module")
def reference() -> dict[str, Verification]:
    """The three official verdicts on `sample_submission/`, replayed through the parser."""
    return {
        task: parse_verifier_output(
            task, (REFERENCE_DIR / f"reference_{task}.json").read_text(), 0.0
        )
        for task in TASKS
    }


@pytest.fixture(scope="module")
def reference_run(reference) -> VerificationRun:
    """Those three verdicts as one finished run."""
    return VerificationRun(results=list(reference.values()))


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
        assert verdict.model_cost == pytest.approx(verdict.score, rel=MODEL_AGREEMENT_REL)
        assert verdict.model_agrees is True

    def test_the_reference_submission_scores_exactly_one(self, reference_run):
        # True by definition -- which makes it a check on the calibration, not the metric.
        costs = reference_run.costs
        assert costs == pytest.approx(DEFAULT_CALIBRATION.reference_costs)
        assert score(costs, DEFAULT_CALIBRATION) == pytest.approx(1.0)

    def test_the_builtin_calibration_names_the_circuits_it_measured(self, reference_run):
        assert DEFAULT_CALIBRATION.reference_sha256 == reference_run.hashes

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
        assert verdict.model_cost is None
        assert verdict.model_agrees is None


class TestDriver:
    """Locating each verifier, pointing it at the right target, and running it."""

    def test_each_task_has_a_script_path(self):
        assert set(VERIFIER_SCRIPTS) == set(TASKS)

    def test_a_missing_verifiers_tree_says_what_to_point_at(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="verifiers"):
            verifier_script(tmp_path, "invsqrt")

    def test_a_verifier_loads_as_a_module_from_its_tree(self, tmp_path, monkeypatch):
        root = fake_verifiers(tmp_path / "verifiers", "dqi", "MAGIC = 41 + 1")
        monkeypatch.setattr(sys, "path", list(sys.path))
        try:
            assert load_verifier(root, "dqi").MAGIC == 42
        finally:
            sys.modules.pop("verify_submission_dqi", None)

    def test_an_unknown_task_is_rejected(self, tmp_path):
        with pytest.raises(KeyError, match="unknown task"):
            verifier_script(tmp_path, "nope")

    def test_a_target_the_audit_did_not_find_is_reported_as_missing(self, tmp_path):
        run = VerificationRun()
        verdicts = list(run.measure(tmp_path, SubmissionAudit.from_root(tmp_path)))
        assert [v.status for v in verdicts] == ["MISSING"] * len(TASKS)
        assert run.results == verdicts

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
        root = fake_verifiers(tmp_path / "verifiers", "invsqrt", ECHO_ARGV)
        verdict = run_verifier(
            root, "invsqrt", tmp_path / "c.kmb", salt="dev-7", seed=1234
        )
        assert verdict.detail["argv"][-4:] == ["--salt", "dev-7", "--seed", "1234"]
        assert (verdict.salt, verdict.seed) == ("dev-7", 1234)

    def test_no_seed_flag_is_passed_when_none_is_asked_for(self, tmp_path):
        root = fake_verifiers(tmp_path / "verifiers", "invsqrt", ECHO_ARGV)
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
        run = VerificationRun(rows_path=path)
        run.record(reference["invsqrt"])
        assert len(path.read_text().strip().splitlines()) == 1
        run.record(reference["dqi"])
        assert len(path.read_text().strip().splitlines()) == 2

    def test_the_rows_read_back_keyed_by_task(self, tmp_path, reference):
        # The provenance record reads this file; writer and reader must agree on its shape.
        path = tmp_path / "rows.jsonl"
        run = VerificationRun(rows_path=path)
        run.record(reference["dqi"])
        read = VerificationRun.read_rows(path)
        assert set(read) == {"dqi"}
        assert read["dqi"]["score"] == pytest.approx(reference["dqi"].score)


class TestMeasuredCalibration:
    """Turning three passing reference verdicts into the calibration -- gate G0's output."""

    def test_passing_verdicts_become_a_measured_calibration(self, reference_run):
        calibration = reference_run.calibration("test")
        assert calibration.measured is True
        assert calibration.reference_costs == pytest.approx(
            DEFAULT_CALIBRATION.reference_costs
        )
        # Without the hashes, calibrating against a compressed candidate would silently make
        # that candidate score 1.0 and every later ratio relative to it.
        assert calibration.reference_sha256 == DEFAULT_CALIBRATION.reference_sha256
        assert calibration.disagreeing_hashes({"dqi": "something-else"}) == {
            "dqi": (DEFAULT_CALIBRATION.reference_sha256["dqi"], "something-else")
        }

    def test_refuses_to_calibrate_from_an_incomplete_run(self, reference):
        failed = Verification("dqi", "FAIL", None, 1.0, detail={"error": "rejected"})
        partial = [reference["invsqrt"], failed, reference["dfthcsos"]]
        with pytest.raises(ValueError, match="must pass to calibrate"):
            VerificationRun(results=partial).calibration("test")
