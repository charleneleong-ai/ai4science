import pytest

from quantum_compilation.cost import (
    MAX_DQI_FAILURE_RATE,
    MIN_E_GAP_HA,
    MIN_LAMBDA,
    PROVISIONAL_CALIBRATION,
    TASKS,
    MetricCalibration,
    dfthcsos_cost,
    dfthcsos_rejections,
    dqi_cost,
    dqi_repetition_factor,
    lambda_eff,
    score,
)

AT_REFERENCE = dict(PROVISIONAL_CALIBRATION.reference_costs)
CLEAN_CHEMISTRY = dict(
    angle_circuit_qubits=[128, 100],
    hamiltonian_error_ha=19e-3,
    e_gap_ha=MIN_E_GAP_HA + 0.1,
    lcu_lambda=MIN_LAMBDA + 1,
)


class TestDqiRepetitionFactor:
    """R(eps) = ln2 / -ln(1/2 + 2 sqrt(eps)), the price of a failing decoder."""

    def test_perfect_decoder_costs_one_repetition(self):
        assert dqi_repetition_factor(0.0) == pytest.approx(1.0)

    def test_grows_with_failure_rate(self):
        factors = [
            dqi_repetition_factor(r) for r in (0.0, 1e-5, 1e-4, 1e-3, 0.01, 0.05)
        ]
        assert factors == sorted(factors)
        assert factors[0] < factors[-1]

    def test_diverges_towards_the_one_sixteenth_pole(self):
        assert dqi_repetition_factor(MAX_DQI_FAILURE_RATE - 1e-12) > 1e9

    @pytest.mark.parametrize("rate", [-1e-9, MAX_DQI_FAILURE_RATE, 0.1, 1.0])
    def test_rejects_rates_where_r_is_not_positive(self, rate):
        with pytest.raises(ValueError, match="failure rate"):
            dqi_repetition_factor(rate)


class TestCostFunctions:
    """The per-task cost formulas, in the terms the Evaluation page states them."""

    def test_dqi_scales_with_the_square_of_the_qubit_count(self):
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
        assert score(improved, PROVISIONAL_CALIBRATION) == pytest.approx(2 ** (1 / 3))

    def test_each_task_carries_equal_weight(self):
        gains = {
            task: score(
                AT_REFERENCE | {task: AT_REFERENCE[task] / 3}, PROVISIONAL_CALIBRATION
            )
            for task in TASKS
        }
        assert gains == pytest.approx({task: 3 ** (1 / 3) for task in TASKS})

    def test_reproduces_the_published_public_frontier_blend(self):
        # The 2026-10-30 briefing (topic 754174) reports per-task ratios of 20.57, 60.59 and
        # 9.93 and a best-per-task geometric mean of 23.13587. Its ratios are quoted to four
        # significant figures, so agreement to ~0.1% is all that can be asked -- but it does
        # confirm the metric is a plain geometric mean of cost ratios, which is the only
        # independent check on the transcription until sample_submission can be measured.
        ratios = {"invsqrt": 20.57, "dqi": 60.59, "dfthcsos": 9.93}
        costs = {task: AT_REFERENCE[task] / r for task, r in ratios.items()}
        assert score(costs, PROVISIONAL_CALIBRATION) == pytest.approx(
            23.13587, rel=1e-3
        )

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
            score(costs, PROVISIONAL_CALIBRATION)

    @pytest.mark.parametrize("cost", [0.0, -1.0])
    def test_rejects_a_non_positive_cost(self, cost):
        with pytest.raises(ValueError, match="must be positive"):
            score(AT_REFERENCE | {"dqi": cost}, PROVISIONAL_CALIBRATION)


class TestChemistryRejections:
    """The dfthcsos constraints that reject a submission outright."""

    def test_a_compliant_submission_has_no_rejections(self):
        assert dfthcsos_rejections(**CLEAN_CHEMISTRY) == []

    @pytest.mark.parametrize(
        ("override", "expected"),
        [
            ({"angle_circuit_qubits": [128, 129]}, "wider than 128 qubits"),
            ({"hamiltonian_error_ha": 21e-3}, "from the reference"),
            ({"lcu_lambda": MIN_LAMBDA}, "Lambda"),
            ({"e_gap_ha": MIN_E_GAP_HA}, "E_gap"),
        ],
    )
    def test_each_constraint_is_reported(self, override, expected):
        reasons = dfthcsos_rejections(**(CLEAN_CHEMISTRY | override))
        assert len(reasons) == 1
        assert expected in reasons[0]

    def test_every_violation_is_listed_not_just_the_first(self):
        reasons = dfthcsos_rejections(
            [200], hamiltonian_error_ha=50e-3, e_gap_ha=1.0, lcu_lambda=1.0
        )
        assert len(reasons) == 4

    def test_an_unsupplied_quantity_is_not_checked_rather_than_assumed_zero(self):
        # Only the widths are known before the verifier runs. A zero E_gap would fail the
        # E_gap > 3.914898 test, so defaulting instead of skipping would invent a rejection.
        assert dfthcsos_rejections(angle_circuit_qubits=[128]) == []


class TestMetricCalibration:
    """The metric's unpublished constants, loaded together with their provenance."""

    def test_the_shipped_calibration_is_marked_unmeasured(self):
        assert PROVISIONAL_CALIBRATION.measured is False
        assert PROVISIONAL_CALIBRATION.rotation_t_cost == {}
        assert "provisional" in PROVISIONAL_CALIBRATION.provenance

    @pytest.mark.parametrize(
        "reference_costs",
        [
            {"invsqrt": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0, "dfthcsos": 1.0, "extra": 1.0},
            {"invsqrt": 1.0, "dqi": 1.0, "typo": 1.0},
        ],
    )
    def test_rejects_a_calibration_that_does_not_cover_every_task(
        self, reference_costs
    ):
        with pytest.raises(ValueError, match="must cover exactly"):
            MetricCalibration(reference_costs, {}, "test", measured=True)

    def test_rejects_a_non_positive_reference_cost(self):
        with pytest.raises(ValueError, match="must be positive"):
            MetricCalibration(AT_REFERENCE | {"dqi": 0.0}, {}, "test", measured=True)

    def test_survives_a_json_round_trip_with_integer_rotation_levels(self, tmp_path):
        # JSON object keys are strings, so the Clifford hierarchy levels have to come back
        # as ints or every pricing lookup silently misses.
        original = MetricCalibration(
            AT_REFERENCE, {4: 9.0, 5: 13.0}, "measured", measured=True
        )
        path = tmp_path / "calibration.json"
        original.to_json(path)
        assert MetricCalibration.from_json(path) == original
