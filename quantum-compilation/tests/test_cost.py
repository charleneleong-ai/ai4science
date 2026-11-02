import json
import math
from functools import partial
from itertools import pairwise
from pathlib import Path
from typing import Any

import pytest
import yaml
from quantum_compilation.cost import (
    BASELINE_FACTORIZATION,
    DEFAULT_CALIBRATION,
    ESCAPE_SCORE,
    AngleFamily,
    AngleFloorSweep,
    and_counting_floor,
    DFTHCSOS_CONSTRAINTS,
    FEMOCO_C0,
    MAX_ANGLE_QUBITS,
    MAX_DQI_FAILURE_RATE,
    MAX_HAMILTONIAN_ERROR_HA,
    MIN_E_GAP_HA,
    MIN_LAMBDA,
    MODEL_AGREEMENT_REL,
    TASKS,
    MetricCalibration,
    clopper_pearson_upper_bound,
    dfthcsos_constraint_margins,
    dfthcsos_cost,
    dfthcsos_rejections,
    dqi_cost,
    dqi_repetition_factor,
    ShapeResources,
    blend_factor,
    cost_cut_for_blend,
    lambda_eff,
    marginal_blend_per_unit,
    lambda_eff_floor,
    lcu_vectors,
    min_lcu_lambda,
    require_task,
    rotation_t_count,
    score,
    word_registers,
    toffoli_equivalents,
)

AT_REFERENCE = dict(DEFAULT_CALIBRATION.reference_costs)

# The plateau submission's leaderboard blend, which several classes price against.
PLATEAU_BLEND = 23.18119


EXPERIMENTS_DIR = Path(__file__).parents[1] / "docs" / "experiments"


@pytest.fixture(scope="module")
def g1() -> dict[str, Any]:
    """Gate G1's committed seed record."""
    return yaml.safe_load((EXPERIMENTS_DIR / "57834335-seeds.yaml").read_text())


def failed_shots(seed: dict[str, Any], shots: int) -> int:
    """Shots the dqi decoder dropped on one seed, from the rate it reported."""
    return round(seed["dqi_detail"]["empirical_error_rate"] * shots)


def reference_report(task: str) -> dict:
    """One shipped verifier's own output, as captured on `sample_submission/`."""
    return json.loads((Path(__file__).parent / "data" / f"reference_{task}.json").read_text())


@pytest.fixture(scope="module")
def dqi_report() -> dict:
    return reference_report("dqi")


@pytest.fixture(scope="module")
def dfthcsos_report() -> dict:
    return reference_report("dfthcsos")


@pytest.fixture(scope="module")
def reference_resources(dfthcsos_report) -> ShapeResources:
    return ShapeResources.from_report(dfthcsos_report)


class TestRotationPricing:
    """`0.56 log2(1/eps) + 5.3` T gates per rotation, the verifier's synthesis model."""

    def test_each_halving_of_the_error_costs_a_fixed_number_of_t_gates(self):
        # Logarithmic, so the price per extra bit of accuracy is a constant 0.56 T -- which
        # also pins the slope's sign, i.e. that a tighter error costs more.
        steps = [rotation_t_count(2.0**-k) for k in (10, 11, 12, 13)]
        assert [b - a for a, b in pairwise(steps)] == pytest.approx([0.56] * 3)

    @pytest.mark.parametrize("eps", [0.0, -1e-9, 1.0, 2.0])
    def test_rejects_an_error_outside_the_unit_interval(self, eps):
        with pytest.raises(ValueError, match="synthesis error"):
            rotation_t_count(eps)


class TestToffoliEquivalents:
    """C = Toffolis + T/4, with rotations priced into T."""

    @pytest.mark.parametrize(("rotation_count", "rotation_eps"), [(0, None), (0.0, 1e-9)])
    def test_four_t_gates_cost_one_toffoli_and_no_rotations_need_no_error(
        self, rotation_count, rotation_eps
    ):
        assert toffoli_equivalents(10, 4, rotation_count, rotation_eps) == pytest.approx(11.0)

    def test_rotations_without_their_error_raise_rather_than_price_at_zero(self):
        # The error is chosen per circuit from the accuracy margin the simulation leaves
        # over, so there is no default that is not a fabrication.
        with pytest.raises(ValueError, match="synthesis error the verifier chose"):
            toffoli_equivalents(10, 4, 37, None)


class TestClopperPearsonBound:
    """The one-sided failure-rate bound the dqi score is computed from."""

    def test_matches_the_bound_the_verifier_computed(self, dqi_report):
        # The reference decoder fails 0 of 102,400 shots and is still bounded near 8.7e-4,
        # which costs it ~19% in repetitions. An observed-fraction model would say zero.
        failures = dqi_report["value_errors"] + dqi_report["phase_errors"]
        assert clopper_pearson_upper_bound(
            failures, dqi_report["verified_shots"], dqi_report["log2_delta"]
        ) == pytest.approx(dqi_report["effective_error_rate"], rel=MODEL_AGREEMENT_REL)

    def test_grows_with_the_failure_count(self):
        bounds = [clopper_pearson_upper_bound(k, 10_000) for k in (0, 1, 5, 50, 500)]
        assert bounds == sorted(bounds)

    def test_shrinks_as_shots_grow(self):
        assert clopper_pearson_upper_bound(0, 10_000) > clopper_pearson_upper_bound(0, 100_000)

    def test_all_shots_failing_bounds_at_one(self):
        assert clopper_pearson_upper_bound(50, 50) == 1.0

    @pytest.mark.parametrize(("failures", "shots"), [(-1, 100), (101, 100), (0, 0), (0, -5)])
    def test_rejects_an_impossible_failure_count(self, failures, shots):
        with pytest.raises(ValueError, match="failures <= shots"):
            clopper_pearson_upper_bound(failures, shots)

    @pytest.mark.parametrize("log2_delta", [1.0, float("nan"), float("-inf")])
    def test_rejects_a_confidence_that_is_not_a_negative_power_of_two(self, log2_delta):
        with pytest.raises(ValueError, match="log2_delta"):
            clopper_pearson_upper_bound(0, 100, log2_delta)


class TestDqiRepetitionFactor:
    """R(eps) = ln2 / -ln(1/2 + 2 sqrt(eps)), the price of a failing decoder."""

    def test_perfect_decoder_costs_one_repetition(self):
        assert dqi_repetition_factor(0.0) == pytest.approx(1.0)

    def test_grows_with_failure_rate(self):
        factors = [dqi_repetition_factor(r) for r in (0.0, 1e-5, 1e-4, 1e-3, 0.01, 0.05)]
        assert factors == sorted(factors)
        assert factors[0] < factors[-1]

    def test_diverges_towards_the_one_sixteenth_pole(self):
        assert dqi_repetition_factor(MAX_DQI_FAILURE_RATE - 1e-12) > 1e9

    @pytest.mark.parametrize("rate", [-1e-9, MAX_DQI_FAILURE_RATE, 0.1, 1.0])
    def test_rejects_rates_where_r_is_not_positive(self, rate):
        with pytest.raises(ValueError, match="failure rate"):
            dqi_repetition_factor(rate)


class TestCostFunctions:
    """The per-task cost formulas, in the terms the verifiers compute them."""

    def test_dqi_scales_with_the_square_of_the_qubit_count(self):
        # Why a qubit is worth twice as much here as on the other two tasks.
        assert dqi_cost(1776, 1e6, 0.0) / dqi_cost(888, 1e6, 0.0) == pytest.approx(4.0)

    def test_lambda_eff_at_the_rejection_boundary(self):
        # An independent anchor: sqrt(3.914898 * (98 - 3.914898)) at the smallest legal inputs.
        assert lambda_eff(MIN_E_GAP_HA, MIN_LAMBDA) == pytest.approx(19.19, abs=0.01)

    def test_lambda_eff_rejects_a_negative_radicand(self):
        with pytest.raises(ValueError, match="not real"):
            lambda_eff(e_gap_ha=100.0, lcu_lambda=1.0)

    def test_dfthcsos_is_linear_in_toffolis(self):
        base = dfthcsos_cost(430, 1e5, 4.0, 50.0)
        assert dfthcsos_cost(430, 2e5, 4.0, 50.0) == pytest.approx(2 * base)


class TestLambdaEffFloor:
    """`D >= -c0`, so `lambda_eff` has a floor -- the bound that closed the refit."""

    def test_the_femoco_floor(self):
        assert lambda_eff_floor(FEMOCO_C0) == pytest.approx(46.786180, abs=1e-6)

    def test_the_derived_bound_puts_D_at_minus_c0(self):
        """What makes the floor a floor: at `min_lcu_lambda`, `D` sits exactly at `-c0`."""
        assert 2 * min_lcu_lambda(FEMOCO_C0) - MIN_E_GAP_HA == pytest.approx(
            -FEMOCO_C0, rel=1e-15
        )

    @pytest.mark.parametrize(
        ("c0", "expected"),
        [
            (FEMOCO_C0, 281.523686),  # the derived bound dominates, 5.7x MIN_LAMBDA
            (-50.0, MIN_LAMBDA),  # a shallow c0 leaves the verifier's rule binding instead
        ],
    )
    def test_the_one_norm_floor_takes_whichever_bound_binds(self, c0, expected):
        assert min_lcu_lambda(c0) == pytest.approx(expected, abs=1e-6)


class TestScore:
    """The geometric mean of the three per-task cost ratios."""

    def test_halving_one_task_moves_the_score_by_the_cube_root(self):
        improved = AT_REFERENCE | {"dqi": AT_REFERENCE["dqi"] / 2}
        assert score(improved, DEFAULT_CALIBRATION) == pytest.approx(2 ** (1 / 3))

    def test_each_task_carries_equal_weight(self):
        gains = {
            task: score(AT_REFERENCE | {task: AT_REFERENCE[task] / 3}, DEFAULT_CALIBRATION)
            for task in TASKS
        }
        assert gains == pytest.approx({task: 3 ** (1 / 3) for task in TASKS})

    def test_reproduces_the_published_public_frontier_blend(self):
        # The 2026-10-30 briefing (topic 754174) reports per-task costs of 100,941,
        # 4.0334e11 and 3.3780e8 and a best-per-task geometric mean of 23.13587. Against the
        # measured reference costs that lands within 0.01%, which confirms both the metric's
        # shape and -- since the briefing's numbers are independent of ours -- the
        # calibration. A third-party report, re-verified by them and not by us.
        costs = {"invsqrt": 100_941, "dqi": 4.0334e11, "dfthcsos": 3.3780e8}
        assert score(costs, DEFAULT_CALIBRATION) == pytest.approx(23.13587, rel=1e-3)

    @pytest.mark.parametrize(
        "costs",
        [
            {"invsqrt": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0, "dfthcsos": 1.0, "extra": 1.0},
        ],
    )
    def test_requires_a_cost_for_every_task(self, costs):
        with pytest.raises(ValueError, match="need a cost for each"):
            score(costs, DEFAULT_CALIBRATION)

    @pytest.mark.parametrize("cost", [0.0, -1.0])
    def test_rejects_a_non_positive_cost(self, cost):
        with pytest.raises(ValueError, match="must be positive"):
            score(AT_REFERENCE | {"dqi": cost}, DEFAULT_CALIBRATION)


class TestBlendPricing:
    """Pricing one task's cost change without the other two, and inverting it."""

    def test_agrees_with_the_full_score(self):
        moved = AT_REFERENCE | {"dqi": AT_REFERENCE["dqi"] / 1.37}
        direct = score(moved, DEFAULT_CALIBRATION)
        assert direct == pytest.approx(
            score(AT_REFERENCE, DEFAULT_CALIBRATION)
            * blend_factor(moved["dqi"], AT_REFERENCE["dqi"])
        )

    @pytest.mark.parametrize(
        ("blend", "target"),
        [(PLATEAU_BLEND, 25.0), (1.0, 2.0), (40.0, 41.0), (25.0, PLATEAU_BLEND)],
    )
    def test_the_cut_it_names_lands_on_the_target(self, blend, target):
        cut = cost_cut_for_blend(blend, target)
        assert blend * blend_factor(1.0 - cut, 1.0) == pytest.approx(target)
        assert (cut < 0) == (target < blend)

    def test_the_per_unit_price_is_the_slope_of_the_factor(self):
        """What `shape` prints per Toffoli, against the finite difference it approximates."""
        total = 16266.0
        step = PLATEAU_BLEND * (blend_factor(total - 1, total) - 1)
        assert marginal_blend_per_unit(PLATEAU_BLEND, total) == pytest.approx(
            step, rel=1e-4
        )

    @pytest.mark.parametrize(
        ("fn", "args"),
        [
            (blend_factor, (0.0, 1.0)),
            (blend_factor, (1.0, -1.0)),
            (cost_cut_for_blend, (0.0, 1.0)),
            (cost_cut_for_blend, (1.0, -1.0)),
            (marginal_blend_per_unit, (0.0, 1.0)),
            (marginal_blend_per_unit, (1.0, -1.0)),
        ],
    )
    def test_rejects_non_positive_inputs(self, fn, args):
        with pytest.raises(ValueError, match="must be positive"):
            fn(*args)


class TestRequireTask:
    """One canonical rejection of an unknown task name, shared by every lookup."""

    @pytest.mark.parametrize("task", ["", "INVSQRT", "not_a_task", "dfthc"])
    def test_rejects_anything_not_in_tasks(self, task):
        with pytest.raises(KeyError, match="unknown task"):
            require_task(task)

    def test_returns_the_task_so_it_can_wrap_a_lookup(self):
        assert require_task("dqi") == "dqi"


class TestShapeResourcesFromReport:
    """The T/Q split against the reference submission's own verifier report.

    Runs on a bare clone: the report is committed, so this is the agreement check the
    closed-form path in `chemistry` cannot make without the 122 MB FCIDUMP.
    """

    @pytest.mark.parametrize("field", ["toffoli", "qubits", "heuristic_toffoli"])
    def test_reproduces_the_reported_total(self, reference_resources, dfthcsos_report, field):
        assert getattr(reference_resources, field) == dfthcsos_report["resources"][field]

    def test_the_angle_share_is_most_of_the_references_toffolis(self, reference_resources):
        """85% here against 51% for the plateau: the share is a property of the shape."""
        assert reference_resources.angle_share == pytest.approx(0.8525, abs=5e-4)

    def test_the_workspace_is_a_max_not_a_sum(self, reference_resources):
        """Both families' temporaries are reported; `Q` counts the larger once."""
        res = reference_resources
        assert res.angle_temporary_qubits < res.heuristic_temporary_qubits
        assert res.qubits < (
            res.heuristic_persistent_qubits
            + res.angle_word_qubits
            + res.heuristic_temporary_qubits
            + res.angle_temporary_qubits
        )

    def test_J_is_recovered_from_the_reported_total(self, reference_resources, dfthcsos_report):
        """`word_registers` inverts `Q_total`; the report states J directly, so it is checkable."""
        res, resources = reference_resources, dfthcsos_report["resources"]
        recovered = word_registers(
            resources["qubits"],
            res.heuristic_persistent_qubits,
            res.heuristic_temporary_qubits,
            resources["brot"],
        )
        assert recovered == dfthcsos_report["batch_size"]
        assert recovered * resources["brot"] == res.angle_word_qubits

    @pytest.mark.parametrize("qubits_total", [90, 731])
    def test_an_undecomposable_total_is_rejected(self, qubits_total):
        with pytest.raises(ValueError, match="word registers"):
            word_registers(qubits_total, 90, 544, 12)


class TestChemistryRejections:
    """The one dfthcsos constraint knowable before the verifier runs: the circuit widths."""

    def test_a_compliant_submission_has_no_rejections(self):
        # The width limit is a cap, so sitting exactly on 128 is legal.
        assert dfthcsos_rejections([MAX_ANGLE_QUBITS, 100]) == []

    def test_every_circuit_over_the_cap_is_counted_in_one_reason(self):
        reasons = dfthcsos_rejections([200, 300, 64])
        assert reasons == [
            f"2 ANGLE circuit(s) wider than {MAX_ANGLE_QUBITS} qubits: [200, 300]"
        ]


class TestConstraintTable:
    """The four hard limits, and the cap/floor distinction that decides a value on the line."""

    def test_the_table_covers_exactly_the_constraints_that_reject(self):
        # The rejection rule and the margins read the same table, so a submission cannot be
        # reported as having headroom against a limit nobody checks.
        assert {c.name for c in DFTHCSOS_CONSTRAINTS} == {
            "hamiltonian_error_ha",
            "angle_qubits",
            "e_gap_ha",
            "lcu_lambda",
        }

    @pytest.mark.parametrize(
        ("name", "legal_at_limit"),
        [
            ("hamiltonian_error_ha", True),
            ("angle_qubits", True),
            ("e_gap_ha", False),
            ("lcu_lambda", False),
        ],
    )
    def test_caps_admit_their_limit_and_the_floors_reject_it(self, name, legal_at_limit):
        # The verifier rejects an error or width *above* its cap but rejects E_gap or Lambda
        # at or below its floor: the whole difference a cap and a floor make.
        constraint = next(c for c in DFTHCSOS_CONSTRAINTS if c.name == name)
        assert constraint.legal_at_limit is legal_at_limit
        assert constraint.violated(constraint.limit) is not legal_at_limit
        inside = math.nextafter(constraint.limit, 0.0 if legal_at_limit else 1e9)
        assert not constraint.violated(inside)


class TestFactorizationShape:
    """The (R, B, C) vocabulary the verifier and the parts manifest share."""

    @pytest.mark.parametrize(
        ("shape", "angle_words"),
        [
            # The row count of each submission's npz `angle_words`, read off the files: the
            # shipped sample_submission at the baseline shape, and public submission 57834335
            # at the shape its author declared in the parts manifest. Both are measurements,
            # so this pins the formula rather than restating it.
            (BASELINE_FACTORIZATION[:2], 3976),
            ((8, 143), 1220),
        ],
    )
    def test_the_vector_count_matches_the_submitted_angle_words(self, shape, angle_words):
        assert lcu_vectors(*shape) == angle_words


class TestChemistryConstraintMargins:
    """Closeness-to-rejection for each constraint, read off a verifier report."""

    def test_sitting_on_a_cap_and_on_a_floor_both_read_as_fully_spent(self):
        # A cap is approached from below and a floor from above, so the fraction has to
        # invert for floors; both directions must give 1.0 at the limit.
        margins = dfthcsos_constraint_margins(
            {
                "metrics": {
                    "error_1": 0.0,
                    "error_2": MAX_HAMILTONIAN_ERROR_HA,
                    "E_gap": MIN_E_GAP_HA,
                    "Lambda": MIN_LAMBDA,
                },
                "resources": {"angles_qubits": MAX_ANGLE_QUBITS},
            }
        )
        assert all(m["fraction_of_limit_used"] == 1.0 for m in margins.values())

    def test_the_binding_accuracy_term_is_the_worse_of_the_two(self):
        margins = dfthcsos_constraint_margins(
            {
                "metrics": {
                    "error_1": 19e-3,
                    "error_2": 3e-3,
                    "E_gap": MIN_E_GAP_HA + 1,
                    "Lambda": MIN_LAMBDA + 1,
                },
                "resources": {"angles_qubits": 64},
            }
        )
        assert margins["hamiltonian_error_ha"]["value"] == 19e-3

    def test_the_reference_submission_has_headroom_on_every_constraint(self):
        margins = dfthcsos_constraint_margins(reference_report("dfthcsos"))
        assert all(0 < m["fraction_of_limit_used"] < 1 for m in margins.values())
        # Lambda is the one constraint nowhere near binding -- the distinction that makes
        # MIN_LAMBDA a Lambda floor rather than a lambda_eff floor.
        assert margins["lcu_lambda"]["fraction_of_limit_used"] < 0.2


class TestMetricCalibration:
    """The metric's one unpublished input per task, loaded with its provenance."""

    def test_the_shipped_calibration_is_measured(self):
        assert DEFAULT_CALIBRATION.measured is True
        assert "sample_submission" in DEFAULT_CALIBRATION.provenance

    @pytest.mark.parametrize(
        "reference_costs",
        [
            {"invsqrt": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0, "dfthcsos": 1.0, "extra": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0, "typo": 1.0},
        ],
    )
    def test_rejects_a_calibration_that_does_not_cover_every_task(self, reference_costs):
        with pytest.raises(ValueError, match="must cover exactly"):
            MetricCalibration(reference_costs, "test", measured=True)

    def test_rejects_a_non_positive_reference_cost(self):
        with pytest.raises(ValueError, match="must be positive"):
            MetricCalibration(AT_REFERENCE | {"dqi": 0.0}, "test", measured=True)

    def test_rejects_hashes_for_a_task_that_does_not_exist(self):
        with pytest.raises(ValueError, match="non-task key"):
            MetricCalibration(AT_REFERENCE, "test", measured=True, reference_sha256={"x": "a"})

    def test_survives_a_json_round_trip(self, tmp_path):
        path = tmp_path / "calibration.json"
        DEFAULT_CALIBRATION.to_json(path)
        assert MetricCalibration.from_json(path) == DEFAULT_CALIBRATION

class TestSeedNoiseRecord:
    """Gate G1's record, re-derived from the cost model that quotes it."""

    def test_the_dqi_bound_and_repetition_factor_follow_from_the_failure_count(self, g1):
        shots = g1["shots_per_run"]
        for seed in g1["seeds"]:
            dqi = seed["dqi_detail"]
            assert clopper_pearson_upper_bound(
                failed_shots(seed, shots), shots
            ) == pytest.approx(dqi["effective_error_rate"], rel=MODEL_AGREEMENT_REL)
            assert dqi_repetition_factor(dqi["effective_error_rate"]) == pytest.approx(
                dqi["repetition_factor"], rel=MODEL_AGREEMENT_REL
            )

    def test_the_whole_blend_spread_comes_from_dqi(self, g1):
        """The claim the README and the artifact both make: chemistry carries no noise."""
        spread = {
            task: max(s["costs"][task] for s in g1["seeds"])
            - min(s["costs"][task] for s in g1["seeds"])
            for task in TASKS
        }
        assert spread["dfthcsos"] == 0.0
        assert spread["invsqrt"] / g1["seeds"][0]["costs"]["invsqrt"] < 1e-6
        assert spread["dqi"] / g1["seeds"][0]["costs"]["dqi"] > 1e-3
        counts = {failed_shots(s, g1["shots_per_run"]) for s in g1["seeds"]}
        assert counts == {0, 1}, "the measured step is 0-or-1 failures; a wider set is a new fact"

    def test_the_blend_moves_only_with_the_failure_count(self, g1):
        """Seeds that dropped no shot must share a blend; the one that dropped one must not."""
        blends: dict[int, set[float]] = {}
        for seed in g1["seeds"]:
            count = failed_shots(seed, g1["shots_per_run"])
            blends.setdefault(count, set()).add(
                round(score(seed["costs"], DEFAULT_CALIBRATION), 4)
            )
        assert set(blends) == {0, 1}, "both a clean and a failed seed are needed to compare"
        assert len(blends[0]) == 1, "zero-failure seeds disagree on the blend"
        assert blends[0] != blends[1], "the failed seed scored the same as the clean ones"


# The plateau submission 57834335. The heuristic split is not in the parts record, so those
# five terms are transcribed -- but the three totals they must reproduce are, and
# `test_the_fixture_still_matches_the_record` is what catches drift in either direction.
def output_bit_floor(input_bits: int, num_addresses: int, output_bits: int) -> int:
    """`and_counting_floor` with its arguments ordered for a sweep over the output width."""
    return and_counting_floor(output_bits, input_bits + output_bits, num_addresses)


PARTS = yaml.safe_load((EXPERIMENTS_DIR / "57834335-parts.yaml").read_text())["parts"]["dfthcsos"]
INCUMBENT = ShapeResources(
    outer_toffoli=369,
    inner_toffoli=5311,
    givens_toffoli=2250,
    heuristic_persistent_qubits=86,
    heuristic_temporary_qubits=304,
    angle_toffoli=float(PARTS["T_angles"]),
    angle_word_qubits=40,
    angle_temporary_qubits=75,
)
SUBMITTED_FAMILY = AngleFamily.from_resources(
    INCUMBENT, int(PARTS["brot"]), lcu_vectors(int(PARTS["R"]), int(PARTS["B"]))
)
SUBMITTED_ANDS = SUBMITTED_FAMILY.ands_for_toffoli(INCUMBENT.angle_toffoli)


# The sweep as the CLI builds it, off the five transcribed heuristic terms, and as the artifact
# builder builds it, off the committed record alone. They must agree to the last digit: the
# page has no verifier on its path and so cannot check its own numbers.
SWEEP = AngleFloorSweep.over(
    INCUMBENT, int(PARTS["brot"]), SUBMITTED_FAMILY.num_addresses, PLATEAU_BLEND
)
RECORD_SWEEP = AngleFloorSweep.over(
    ShapeResources.from_parts_record(PARTS),
    int(PARTS["brot"]),
    SUBMITTED_FAMILY.num_addresses,
    PLATEAU_BLEND,
)


class TestAngleLoaderFloor:
    """The counting floor on the loaders, and the batch size it turns into a lever."""

    def test_the_fixture_still_matches_the_record(self):
        """The five transcribed heuristic terms, checked against the totals that are recorded."""
        assert INCUMBENT.toffoli == PARTS["T_total"]
        assert INCUMBENT.qubits == PARTS["Q_total"]
        assert SUBMITTED_FAMILY.num_circuits == 15
        assert SUBMITTED_FAMILY.address_bits == 11

    def test_the_recorded_batch_size_agrees_with_what_the_total_implies(self):
        """`J` is ours, written by the generator, so the third party's `Q_total` audits it."""
        assert int(PARTS["J"]) == SUBMITTED_FAMILY.batch_size
        assert int(PARTS["J"]) == word_registers(
            int(PARTS["Q_total"]),
            INCUMBENT.heuristic_persistent_qubits,
            INCUMBENT.heuristic_temporary_qubits,
            int(PARTS["brot"]),
        )

    def test_the_submitted_loaders_sit_above_the_floor_but_well_under_twice_it(self):
        """What makes this a floor worth quoting: the loaders are already close to it."""
        assert SUBMITTED_FAMILY.and_floor == 157
        assert 1.0 < SWEEP.rate < 2.0

    def test_granting_the_word_registers_as_inputs_lowers_the_floor(self):
        """They arrive holding a function of the address, so a bound must allow their use."""
        assert SUBMITTED_FAMILY.and_floor < SUBMITTED_FAMILY.address_only_floor == 191

    @pytest.mark.parametrize(
        ("floor_of", "sizes"),
        [
            (partial(and_counting_floor, 40, 51), range(200, 2001, 200)),
            (partial(output_bit_floor, 11, 1220), range(8, 129, 8)),
        ],
        ids=["addresses", "output_bits"],
    )
    def test_the_floor_grows_with_the_table_it_has_to_encode(self, floor_of, sizes):
        assert all(a < b for a, b in pairwise([floor_of(n) for n in sizes]))

    @pytest.mark.parametrize("bits", [(0, 11, 1220), (40, 0, 1220), (40, 11, 0), (-40, 11, 1220)])
    def test_a_degenerate_shape_is_rejected_rather_than_returning_a_floor(self, bits):
        with pytest.raises(ValueError, match="positive"):
            and_counting_floor(*bits)

    def test_the_width_cap_is_what_bounds_the_batch_size(self):
        """J=14 is the widest legal batching at brot=8: 11 + 112 fits 128, J=15 does not."""
        brot, addresses = SUBMITTED_FAMILY.brot, SUBMITTED_FAMILY.num_addresses
        assert AngleFamily(14, brot, addresses).feasible
        assert not AngleFamily(15, brot, addresses).feasible
        assert [row.batch_size for row in SWEEP.rows] == list(range(1, 15))

    def test_re_batching_the_submission_into_its_own_family_changes_nothing(self):
        """The round trip that keeps the projection honest, exact on every resource term."""
        same = INCUMBENT.with_angle_family(
            SUBMITTED_FAMILY, SUBMITTED_ANDS, INCUMBENT.angle_temporary_qubits
        )
        assert same == INCUMBENT
        assert SWEEP.at(SUBMITTED_FAMILY.batch_size).at_rate == pytest.approx(PLATEAU_BLEND)

    def test_the_floor_permits_the_escape_score_once_the_word_inputs_are_granted(self):
        """The verdict r3 turns on: the loaders are not closed off from the escape band."""
        assert SWEEP.permits_target
        assert SWEEP.best_floor.batch_size == 11
        assert SWEEP.best_floor.at_floor == pytest.approx(25.684, abs=5e-4)
        assert SWEEP.best_floor.at_floor > ESCAPE_SCORE

    def test_re_batching_alone_is_worth_a_point_at_the_submitted_and_rate(self):
        """J is a free parameter, so this gain needs no new circuit cleverness."""
        best = SWEEP.best_rate
        assert best.batch_size == 13
        assert best.at_rate - PLATEAU_BLEND == pytest.approx(1.0405, abs=5e-4)
        assert best.at_rate < ESCAPE_SCORE, "re-batching alone is not the escape"
        assert best.target_over_floor == pytest.approx(1.326, abs=5e-4)

    def test_the_record_alone_prices_the_sweep_exactly_as_the_split_does(self):
        """What licenses the artifact page: it has no verifier, so the lumping must be free.

        Folding the heuristic workspace into the persistent count is faithful only while that
        workspace is the larger side of the `max` in `qubits`, which the ANGLE width cap
        guarantees -- so the second assertion is the reason the first one holds.
        """
        for lumped, split in zip(RECORD_SWEEP.rows, SWEEP.rows, strict=True):
            assert lumped.batch_size == split.batch_size
            assert lumped.at_floor == pytest.approx(split.at_floor, abs=1e-12)
            assert lumped.at_rate == pytest.approx(split.at_rate, abs=1e-12)
            assert lumped.ands_for_target == pytest.approx(split.ands_for_target, abs=1e-9)
        assert max(row.family.free_workspace for row in SWEEP.rows) < (
            INCUMBENT.heuristic_temporary_qubits
        )
