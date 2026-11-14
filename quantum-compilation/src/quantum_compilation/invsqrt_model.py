"""A bit-exact model of the seed-plus-series invsqrt design, and what each variant costs.

The decoded circuit (docs/invsqrt-diagnosis.md) kicks a short seed `y0 ~ 1/sqrt(v)` from a table,
then a correction. This models that family as

    eps = 1 - y0^2 v,    y = y0 (1 + eps/2 + 3 eps^2/8 + 5 eps^3/16 ...)

in fixed point with `frac_bits` fractional bits, every product floored, so its phase can be checked
against the verifier's target on every input. Costs count partial products of truncated multipliers
and adder widths; they are estimates, calibrated as a whole against the measured circuit.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .cost import rotation_eps, toffoli_equivalents
from .invsqrt import MANTISSA_BITS

GUARD_BITS = 2  # kept below the cut by truncated multipliers, so their error stays under an ulp
CLIFFORD_KICK_BITS = 3  # phases of 2^-1..2^-3 turns are S/T-level, not fine rotations
TAYLOR = (1 / 2, 3 / 8, 5 / 16)  # (1 - eps)^(-1/2) = 1 + eps/2 + 3 eps^2/8 + 5 eps^3/16 + ...
CONSTANT_BITS = 16  # fractional bits of the series coefficients
FIT_NODES = 64
REGISTER_QUBITS = MANTISSA_BITS + 2  # the input, the exponent parity, one lookup scratch bit


def ceil_shift(values: np.ndarray, shift: np.ndarray | int) -> np.ndarray:
    return -((-values) >> shift)


def bit_width(values: np.ndarray) -> int:
    """Bits for the largest magnitude, plus a sign bit if both signs occur."""
    signed = bool((values < 0).any() and (values > 0).any())
    return int(np.abs(values).max()).bit_length() + signed


def partial_products(a: tuple[int, int], b: tuple[int, int], cut: int) -> int:
    """Partial products of operands spanning bit positions `[lo, hi]` kept at or above `cut`."""
    return sum(max(0, b[1] - max(b[0], cut - i) + 1) for i in range(a[0], a[1] + 1))


def square_products(a: tuple[int, int], cut: int) -> int:
    """A squarer needs each cross term once: half the off-diagonal products plus the diagonal."""
    diagonal = sum(2 * i >= cut for i in range(a[0], a[1] + 1))
    return (partial_products(a, a, cut) - diagonal) // 2 + diagonal


def signed_digits(value: int) -> int:
    """Nonzero digits of `value` in non-adjacent form: the adds a constant multiply costs, plus one."""
    digits = 0
    while value:
        if value & 1:
            value -= 2 - (value & 3)  # NAF digit +1 or -1
            digits += 1
        value >>= 1
    return digits


def fixed_point_span(frac_bits: int, width: int) -> tuple[int, int]:
    """Bit positions `[lo, hi]` of a `width`-bit value stored at 2^-frac_bits."""
    return -frac_bits, width - 1 - frac_bits


def floored(phase: np.ndarray, kick_floor: int) -> np.ndarray:
    """The phase with every bit below 2^-kick_floor turns dropped (units of 2^-64 turns)."""
    return phase & ~np.uint64(2 ** (64 - kick_floor) - 1)


@dataclass(frozen=True)
class SeriesDesign:
    """One variant: table size, seed width, series order, precision, and what the table holds."""

    address_bits: int  # mantissa bits after the leading one that index the table (plus parity)
    seed_bits: int  # fractional bits of y0
    terms: int  # series terms after the 1
    frac_bits: int  # fractional bits of eps and the correction
    kick_floor: int  # smallest kicked phase is 2^-kick_floor turns
    store_offset: bool  # table holds 1 - y0^2 v_segment; else eps comes from the full mantissa
    store_square: bool  # table holds y0^2 2^parity; else it is squared on the fly
    tuned: bool = False  # series coefficients fitted to the eps range; else Taylor's

    @property
    def low_bits(self) -> int:
        return MANTISSA_BITS - 1 - self.address_bits

    @cached_property
    def keys(self) -> tuple[np.ndarray, np.ndarray]:
        """Each table entry's segment and exponent parity, indexed by `segment | parity << A`."""
        key = np.arange(2 ** (self.address_bits + 1))
        return key & (2**self.address_bits - 1), key >> self.address_bits

    @cached_property
    def seeds(self) -> np.ndarray:
        """`round(2^S / sqrt(v))` at each segment's midpoint."""
        segment, parity = self.keys
        v = (1 + (segment + 0.5) / 2**self.address_bits) * 2.0**parity
        return np.round(2.0**self.seed_bits / np.sqrt(v)).astype(np.int64)

    @cached_property
    def coefficients(self) -> np.ndarray:
        """Series coefficients at 2^-CONSTANT_BITS: Taylor's, or a least-squares fit at Chebyshev
        nodes over the eps every segment can reach, which absorbs the first dropped term."""
        if not self.tuned:
            return np.round(np.array(TAYLOR[: self.terms]) * 2**CONSTANT_BITS).astype(np.int64)
        segment, parity = self.keys
        square = (self.seeds / 2.0**self.seed_bits) ** 2 * 2.0**parity
        lo = 1 - square * (1 + (segment + 1) / 2**self.address_bits)
        hi = 1 - square * (1 + segment / 2**self.address_bits)
        nodes = np.cos(np.pi * (np.arange(FIT_NODES) + 0.5) / FIT_NODES)
        eps = lo.min() + (hi.max() - lo.min()) * (nodes + 1) / 2
        basis = eps[:, None] ** np.arange(1, self.terms + 1)
        fit, *_ = np.linalg.lstsq(basis, (1 - eps) ** -0.5 - 1, rcond=None)
        return np.round(fit * 2**CONSTANT_BITS).astype(np.int64)

    def evaluate(self, mantissa: np.ndarray, exponent: np.ndarray) -> tuple[Evaluation, np.ndarray]:
        """The value widths, and each input's phase in 2^-64 turns before the kick floor."""
        parity, half = exponent % 2, exponent // 2
        s, f = self.seed_bits, self.frac_bits
        segment = (mantissa >> self.low_bits) & (2**self.address_bits - 1)
        y0 = self.seeds[segment | parity << self.address_bits]
        square = y0 * y0 << parity  # y0^2 2^p at 2^-2S
        shift = 2 * s + MANTISSA_BITS - 1 - f  # v = M 2^p / 2^21; products land at 2^-F
        offset_width = 0
        if self.store_offset:
            base = (1 << (MANTISSA_BITS - 1)) + (segment << self.low_bits)
            offset = (1 << f) - ceil_shift(square * base, shift)
            eps = offset - (square * (mantissa & (2**self.low_bits - 1)) >> shift)
            offset_width = bit_width(offset)
        else:
            eps = (1 << f) - (square * mantissa >> shift)
        power, z, powers = eps, np.zeros_like(eps), []
        for k, coefficient in enumerate(self.coefficients):
            if k:
                power = power * eps >> f
            powers.append(bit_width(power))
            z += power * coefficient >> CONSTANT_BITS
        correction = y0 * z >> s
        y = (y0 << (f - s)) + correction
        phase = np.left_shift(y.astype(np.uint64), (64 - f - half).astype(np.uint64))
        widths = Widths(offset_width, bit_width(eps), powers, bit_width(z), bit_width(correction))
        return Evaluation(self, widths, int(half.max())), phase


@dataclass(frozen=True)
class Widths:
    """Bits each fixed-point value needs over the inputs, sign included where both signs occur."""

    offset: int
    eps: int
    powers: list[int]
    z: int
    correction: int


@dataclass(frozen=True)
class CostModel:
    """The two constants the circuit does not pin down; the search reports how much they matter."""

    toffolis_per_partial_product: float = 2.0  # an AND to form it plus an adder Toffoli
    lookup_toffolis_per_entry: float = 129 / 256  # 6c6a734d's lookup: 129 CCX for 256 entries


@dataclass(frozen=True)
class Evaluation:
    """What costing a design needs from running it: its value widths."""

    design: SeriesDesign
    widths: Widths
    max_half_exponent: int

    def at_floor(self, kick_floor: int) -> Evaluation:
        """The same arithmetic kicked down to another floor: only the rotation count changes."""
        return dataclasses.replace(self, design=dataclasses.replace(self.design, kick_floor=kick_floor))

    @property
    def table_bits(self) -> int:
        d = self.design
        return d.seed_bits + (2 * d.seed_bits + 1) * d.store_square + self.widths.offset

    @property
    def rotations(self) -> int:
        """Fine kicks: the seed shifted by every half exponent, then the correction down to the floor."""
        d = self.design
        seed = d.seed_bits + self.max_half_exponent - CLIFFORD_KICK_BITS
        top = self.widths.correction - 1 - d.frac_bits  # highest correction bit at 2^top
        return seed + d.kick_floor + top + 1

    def stage_toffolis(self, model: CostModel) -> dict[str, float]:
        d, w = self.design, self.widths
        f, s = d.frac_bits, d.seed_bits
        cut = -f - GUARD_BITS
        square = (-2 * s, 0)  # y0^2 2^p
        mantissa = (-(MANTISSA_BITS - 1), d.low_bits - MANTISSA_BITS if d.store_offset else 0)
        eps = fixed_point_span(f, w.eps)
        products = {
            "square seed": 0 if d.store_square else square_products((-s, -1), -2 * s),
            "eps product": partial_products(square, mantissa, cut),
            "eps powers": sum(
                square_products(eps, cut)
                if k == 1
                else partial_products(fixed_point_span(f, w.powers[k - 1]), eps, cut)
                for k in range(1, d.terms)
            ),
            "correction product": partial_products((-s, -1), fixed_point_span(f, w.z), cut),
        }
        adds = sum(signed_digits(int(c)) for c in d.coefficients) - 1 + d.store_offset
        return {
            "lookup": model.lookup_toffolis_per_entry * 2 ** (d.address_bits + 1),
            **{name: n * model.toffolis_per_partial_product for name, n in products.items()},
            "adds": adds * w.z,
            "kicks": self.rotations,  # one CCX loads each kicked bit
        }

    @property
    def work_qubits(self) -> int:
        """Registers live through the arithmetic: seed, eps and the running sum (each power is
        computed, consumed and uncomputed in turn)."""
        return REGISTER_QUBITS + self.design.seed_bits + self.widths.eps + self.widths.z


REFERENCE = SeriesDesign(7, 8, 3, 27, 36, True, False)  # the model's nearest feasible twin of 6c6a734d


@dataclass(frozen=True)
class Calibration:
    """Offsets fitted to one measured circuit (which `REFERENCE` must describe), applied to all."""

    model: CostModel
    toffolis: float  # measured Toffolis the model does not cover (shifter, tail uncompute, ...)
    work_qubits: int

    @classmethod
    def against(
        cls, reference: Evaluation, model: CostModel, toffolis: float, qubits: int
    ) -> Calibration:
        modelled = sum(reference.stage_toffolis(model).values())
        return cls(model, toffolis - modelled, qubits - reference.work_qubits)

    def qubits(self, evaluation: Evaluation) -> int:
        """The lookup peak (register + table, which the measured circuit matches exactly) or the
        calibrated arithmetic peak, whichever is larger."""
        return max(REGISTER_QUBITS + evaluation.table_bits, evaluation.work_qubits + self.work_qubits)

    def cost(self, evaluation: Evaluation, margin: int) -> dict[str, float]:
        n_rot = evaluation.rotations
        eps = rotation_eps(margin, n_rot)
        toffolis = self.toffolis + sum(evaluation.stage_toffolis(self.model).values())
        c = toffoli_equivalents(toffolis, 0, n_rot, eps)
        q = self.qubits(evaluation)
        return {"qubits": q, "toffolis": toffolis, "C": c, "cost": q * c, "eps": eps, "rotations": n_rot}
