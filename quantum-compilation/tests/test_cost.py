import json
import math
from itertools import pairwise
from pathlib import Path

import pytest
from quantum_compilation.cost import (
    BASELINE_FACTORIZATION,
    DEFAULT_CALIBRATION,
    DFTHCSOS_CONSTRAINTS,
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
    lambda_eff,
    lcu_vectors,
    require_task,
    rotation_t_count,
    score,
    toffoli_equivalents,
)

AT_REFERENCE = dict(DEFAULT_CALIBRATION.reference_costs)


@pytest.fixture(scope="module")
def dqi_report() -> dict:
    """The official dqi verifier's own output, as captured on `sample_submission/`."""
    return json.loads((Path(__file__).parent / "data" / "reference_dqi.json").read_text())


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


class TestRequireTask:
    """One canonical rejection of an unknown task name, shared by every lookup."""

    @pytest.mark.parametrize("task", ["", "INVSQRT", "not_a_task", "dfthc"])
    def test_rejects_anything_not_in_tasks(self, task):
        with pytest.raises(KeyError, match="unknown task"):
            require_task(task)

    def test_returns_the_task_so_it_can_wrap_a_lookup(self):
        assert require_task("dqi") == "dqi"


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
        detail = json.loads(
            (Path(__file__).parent / "data" / "reference_dfthcsos.json").read_text()
        )
        margins = dfthcsos_constraint_margins(detail)
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
