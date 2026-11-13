"""Inputs and number formats shared by the invsqrt diagnostic scripts."""

from __future__ import annotations

import numpy as np

MANTISSA_BITS = 22  # x normalised by a power of two to [2^21, 2^22)
DENSE_INPUTS = 4096  # every x up to this is always sampled


def sample_inputs(max_x: int, count: int, seed: int) -> np.ndarray:
    """Every x up to 4096, both neighbours of each power of two, and the rest uniform.

    Not the verifier's `sample_shots` sweep: the normalising shift is where circuits break, so
    its edges matter more than matching the verifier's order.
    """
    edges = [v for k in range(1, max_x.bit_length()) for v in (2**k - 1, 2**k, 2**k + 1)]
    fixed = np.unique(np.r_[np.arange(1, DENSE_INPUTS + 1), edges, max_x].clip(1, max_x))
    rest = np.random.default_rng(seed).integers(1, max_x + 1, max(count - len(fixed), 0))
    return np.r_[fixed, rest].astype(np.uint64)


def normalised(xs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`x * 2^s` in [2^21, 2^22), rounded down, and the exponent `floor(log2 x)`."""
    exponent = np.floor(np.log2(xs.astype(np.float64))).astype(np.int64)
    mantissa = np.ldexp(xs.astype(np.float64), MANTISSA_BITS - 1 - exponent)  # exact below 2^53
    return np.floor(mantissa).astype(np.int64), exponent
