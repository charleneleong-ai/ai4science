"""Input sampling and the power-of-two normalisation the invsqrt scripts share."""

import numpy as np
import pytest
from quantum_compilation.invsqrt import DENSE_INPUTS, MANTISSA_BITS, normalised, sample_inputs

MAX_X = 3_133_452


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
