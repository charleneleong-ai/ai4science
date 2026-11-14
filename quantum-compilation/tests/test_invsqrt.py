"""Input sampling, normalisation, and the seed-plus-series design model for invsqrt."""

import dataclasses
import math

import numpy as np
import pytest
from quantum_compilation.cost import PHASE_UNITS_PER_TURN
from quantum_compilation.invsqrt import (
    DENSE_INPUTS,
    MANTISSA_BITS,
    normalised,
    phase_margin,
    sample_inputs,
)
from quantum_compilation.invsqrt_model import (
    REFERENCE,
    Calibration,
    CostModel,
    SeriesDesign,
    floored,
    partial_products,
    signed_digits,
    square_products,
)

MAX_X = 3_133_452


def reference_values(xs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The verifier's `target` and `tol` for beta = 1, precision 24, recomputed independently so
    these tests run without the competition download."""
    rs = [math.isqrt((1 << 130) // x) for x in xs.tolist()]
    target = np.array([((r + 1) >> 1) % PHASE_UNITS_PER_TURN for r in rs], dtype=np.uint64)
    return target, np.array([r >> 25 for r in rs], dtype=np.uint64)


@pytest.fixture(scope="module")
def sample() -> tuple[tuple[np.ndarray, np.ndarray], np.ndarray, np.ndarray]:
    xs = sample_inputs(MAX_X, 20_000, seed=0)
    return normalised(xs), *reference_values(xs)


def margin(design: SeriesDesign, sample) -> int:
    inputs, target, tol = sample
    _, phase = design.evaluate(*inputs)
    return phase_margin(floored(phase, design.kick_floor), target, tol)


class TestSampleInputs:
    def test_covers_the_dense_range_every_power_of_two_edge_and_the_maximum(self):
        xs = set(sample_inputs(MAX_X, 20_000, seed=0).tolist())
        edges = {v for k in range(1, MAX_X.bit_length()) for v in (2**k - 1, 2**k, 2**k + 1)}
        assert set(range(1, DENSE_INPUTS + 1)) | {v for v in edges if v <= MAX_X} | {MAX_X} <= xs

    def test_stays_in_range_and_fills_up_to_the_count(self):
        xs = sample_inputs(MAX_X, 20_000, seed=0)
        assert len(xs) == 20_000 and xs.min() >= 1 and xs.max() <= MAX_X


class TestNormalised:
    @pytest.mark.parametrize("x", [1, 2, 3, 4096, 4097, 2**21 - 1, 2**21, MAX_X])
    def test_mantissa_is_x_shifted_into_the_top_binade(self, x):
        mantissa, exponent = normalised(np.array([x], dtype=np.uint64))
        assert exponent[0] == x.bit_length() - 1
        assert mantissa[0] == x << (MANTISSA_BITS - x.bit_length())


class TestCounting:
    @pytest.mark.parametrize(
        ("value", "digits"), [(0, 0), (1, 1), (3, 2), (5, 2), (7, 2), (15, 2), (11, 3)]
    )
    def test_signed_digits_follow_the_non_adjacent_form(self, value, digits):
        assert signed_digits(value) == digits

    def test_a_cut_drops_the_low_partial_products(self):
        assert partial_products((0, 1), (0, 1), cut=-10) == 4
        assert partial_products((0, 1), (0, 1), cut=1) == 3

    def test_a_square_counts_each_cross_term_once(self):
        assert square_products((0, 1), cut=-10) == 3  # a0^2, a1^2, a0 a1


class TestSeriesDesign:
    @pytest.mark.parametrize(("terms", "passes"), [(1, False), (2, False), (3, True)])
    def test_an_8_bit_seed_needs_three_terms_to_meet_the_tolerance(self, sample, terms, passes):
        design = SeriesDesign(7, 8, terms, 28, 38, store_offset=True, store_square=False)
        assert (margin(design, sample) > 0) is passes

    def test_a_10_bit_seed_meets_it_with_two(self, sample):
        assert margin(SeriesDesign(8, 10, 2, 26, 36, store_offset=True, store_square=False), sample) > 0

    def test_a_coarser_kick_floor_only_costs_margin(self, sample):
        coarser = dataclasses.replace(REFERENCE, kick_floor=35)
        assert margin(coarser, sample) < margin(REFERENCE, sample)

    def test_taylor_coefficients_are_exact(self):
        assert REFERENCE.coefficients.tolist() == [2**15, 3 * 2**13, 5 * 2**12]

    def test_the_default_per_product_constant_is_fitted_to_the_reference(self, sample):
        evaluation, _ = REFERENCE.evaluate(*sample[0])
        fitted = CostModel.fitted(evaluation).toffolis_per_partial_product
        assert fitted == pytest.approx(CostModel().toffolis_per_partial_product, abs=0.005)

    def test_calibration_reproduces_the_measured_circuit(self, sample):
        evaluation, _ = REFERENCE.evaluate(*sample[0])
        calibration = Calibration.against(evaluation, CostModel(), toffolis=1218.83, qubits=62)
        cost = calibration.cost(evaluation, margin(REFERENCE, sample))
        assert cost["toffolis"] == pytest.approx(1218.83) and cost["qubits"] == 62
