"""Inputs, number formats and the verifier's phase check, shared by the invsqrt scripts."""

from __future__ import annotations

import functools
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np

from .verify import load_verifier

MANTISSA_BITS = 22  # x normalised by a power of two to [2^21, 2^22)
DENSE_INPUTS = 4096  # every x up to this is always sampled
# Stage starts (instruction index) decoded for 6c6a734d only (docs/invsqrt-diagnosis.md); any other
# circuit needs its own.
STAGES = {
    "shifter": 0,
    "lookup": 1663,
    "lookup cleanup": 7494,
    "pass 1": 8590,
    "arithmetic": 8970,
    "pass 2": 16430,
    "uncompute": 16750,
    "clifford": 25500,
    "tail": 33000,
}


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


def phase_margin(phase: np.ndarray, target: np.ndarray, tolerance: np.ndarray) -> int:
    """The verifier's `min(tol - |error|)` in 2^-64 turns; negative means some input fails."""
    error = (phase - target).view(np.int64)
    return int((tolerance.astype(np.int64) - np.abs(error)).min())


@dataclass(frozen=True)
class ReferenceValues:
    """Every x with the verifier's target phase and tolerance, in units of 2^-64 turns."""

    xs: np.ndarray
    target: np.ndarray
    tol: np.ndarray

    @staticmethod
    @functools.cache
    def load(cache: Path, verifiers: Path) -> ReferenceValues:
        """Tabulated by the verifier once (seconds), then read from `cache`."""
        if not cache.exists():
            vf = load_verifier(verifiers, "invsqrt")
            spec = vf.ProblemSpecInvSqrt.for_n_p(10)
            xs = np.arange(1, spec.max_x + 1, dtype=np.uint64)
            target, tol = vf.reference_values(spec, xs)
            cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cache, xs=xs, target=target, tol=tol)
        data = np.load(cache)
        return ReferenceValues(data["xs"], data["target"], data["tol"])

    @cached_property
    def mantissa_exponent(self) -> tuple[np.ndarray, np.ndarray]:
        return normalised(self.xs)

    def margin(self, phase: np.ndarray) -> int:
        return phase_margin(phase, self.target, self.tol)
