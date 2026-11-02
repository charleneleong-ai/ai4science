"""The competition's three cost functions and the score that combines them.

    Score = (prod_b  cost_ref_b / cost_b) ** (1/3)

    invsqrt   cost = Q * C                with C = Toffolis + T/4
    dqi       cost = R(eps_plus) * Q^2 * T
    dfthcsos  cost = lambda_eff * T * Q   with lambda_eff = sqrt(E_gap (2 Lambda - E_gap))

Transcribed from the **shipped verifiers**, not the Evaluation page, which states the
formulas but not their constants. Three things it omits, each of which changes the number
(the README has the detail and the measured consequences):

1. A rotation costs `0.56 log2(1/eps) + 5.3` T, at a single circuit-wide `eps` the verifier
   picks from the accuracy budget the simulated error leaves over -- not per Clifford level.
2. `T` is the mean Toffoli-like gates executed per shot, not the static gate count.
3. `eps_plus` is a Clopper-Pearson upper bound at `delta = 2^-128`, not the observed rate.

Each of the three is why nothing here can price a circuit without simulating it, which is
what `verify` is for. This module owns both the formulas and the names of the verifier
fields they read (`cost_from_report`), so the model lives in one place.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Literal

from scipy.stats import beta

TASKS = ("invsqrt", "dqi", "dfthcsos")

# C = Toffolis + T/4, from the invsqrt verifier's cost model.
T_PER_TOFFOLI_EQUIVALENT = 4

# Mixed fallback synthesis, arXiv:2203.10064, as priced by the invsqrt verifier.
ROTATION_T_SLOPE = 0.56
ROTATION_T_OFFSET = 5.3

# The dqi verifier's confidence for the Clopper-Pearson bound: delta = 2^-128.
DEFAULT_LOG2_DELTA = -128.0

# How closely `cost_from_report` must reproduce a verifier's own score to count as agreeing.
# One constant for the CLI's per-task tick and the tests, so they cannot drift apart.
MODEL_AGREEMENT_REL = 1e-12

# Hard constraints from the dfthcsos verifier. A submission that trips one is rejected
# outright, not merely scored badly. Both floors are strict -- the verifier rejects on
# `value <= floor` -- and the width cap applies to the widest ANGLE circuit in the family.
MAX_ANGLE_QUBITS = 128
MAX_HAMILTONIAN_ERROR_HA = 20e-3
# MIN_LAMBDA floors the one-norm Lambda, NOT lambda_eff -- the `lcu_lambda` row of
# DFTHCSOS_CONSTRAINTS below reads `("metrics", "Lambda")`, and nothing here ever compares
# lambda_eff to it. The reference's lambda_eff of 49.30
# sits 0.6% above this number by coincidence, which made the two readings indistinguishable;
# public submission 57834335 settles it, passing with lambda_eff 48.293 (below 49.0) and
# Lambda 299.82. Reading it as a lambda_eff floor would invent a cliff that is 530% away.
MIN_LAMBDA = 49.0
MIN_E_GAP_HA = 3.914898

# c0 = E_upper - ecore + n * beta_1, the additive constant of the SOS gap. Every term is fixed
# by the FeMoCo reference and the frozen symmetry shift -- `UPPER_ENERGY = -22140.27335` and
# `ecore` are constants of the verifier's `scoreinfra/reference.py`, and `beta_1` comes from
# its `load_fixed_shifts()`, which overwrites whatever shift a submission carries. It is what
# floors the chemistry score: see `min_lcu_lambda`.
FEMOCO_C0 = -559.132473148504

# The orbital count of the FeMoCo instance every dfthcsos submission factorizes, and the
# baseline factorization's shape (R, B, C) -- the rank, the basis count and the ensemble rank.
# Both come from the shipped verifier rather than from its JSON report, which is why reading
# the report alone cannot recover them: `scoreinfra/reference.py` declares
# `BASELINE_SHAPE = (65, 60, 60)` and documents the submitted arrays as `U` (R, B, N) and
# `W` (R, B, C). The sample submission's npz agrees, carrying `U` of (65, 60, 76).
NUM_ORBITALS = 76
BASELINE_FACTORIZATION = (65, 60, 60)

# Angle words per LCU vector, and so the number of ANGLE circuits a batch size J must cover:
# the Givens network that rotates N orbitals takes N - 1 angles.
NUM_ANGLE_WORDS = NUM_ORBITALS - 1

# The rotation precisions worth costing. No verifier rule bounds `brot`; these are wide enough
# to show both walls that do -- below, the angle grid coarsens past MAX_HAMILTONIAN_ERROR_HA,
# and above, the inner preparation's QROAM jumps to a branch that doubles the shared workspace.
MIN_BROT = 4
MAX_BROT = 16

# The blend that clears the plateau: the lowest score among the teams above the 0.6-wide band
# half the field sits in, off docs/experiments/leaderboard-2026-10-31.json. A result number, so
# every consumer reads it from here rather than typing 25.0.
ESCAPE_SCORE = 25.0

# R(eps) is positive only for eps < 1/16: at eps = 1/16 the denominator -ln(1/2 + 2 sqrt(eps))
# hits zero, so a decoder failing more than 6.25% of shots has unbounded cost. The dqi
# verifier rejects eps_plus >= 1/16 outright.
MAX_DQI_FAILURE_RATE = 1 / 16


def require_task(task: str) -> str:
    """The one place an unknown task name is rejected, so the message is the same everywhere."""
    if task not in TASKS:
        raise KeyError(f"unknown task {task!r}; expected one of {TASKS}")
    return task


@dataclass(frozen=True)
class MetricCalibration:
    """Each task's reference cost -- the one metric input not derivable from a circuit.

    The reference submission scores exactly 1.0 by definition, so these come from running
    the shipped verifiers over `sample_submission/`. `reference_sha256` records which
    circuits produced them: a cost is a simulated mean at a seed drawn from the candidate's
    own bytes, so a calibration is a measurement of specific files, not a constant, and
    pointing `--write-calibration` at the wrong submission would otherwise silently make
    *that* submission score 1.0.
    """

    reference_costs: dict[str, float]
    provenance: str
    measured: bool
    reference_sha256: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if set(self.reference_costs) != set(TASKS):
            raise ValueError(
                f"reference costs must cover exactly {TASKS}, "
                f"got {sorted(self.reference_costs)}"
            )
        if bad := {t: c for t, c in self.reference_costs.items() if c <= 0}:
            raise ValueError(f"reference costs must be positive, got {bad}")
        if unknown := sorted(set(self.reference_sha256) - set(TASKS)):
            raise ValueError(f"reference_sha256 has non-task key(s) {unknown}")

    def disagreeing_hashes(self, hashes: Mapping[str, str]) -> dict[str, tuple[str, str]]:
        """Tasks whose circuits differ from the ones this calibration was measured on."""
        return {
            task: (self.reference_sha256[task], hashes[task])
            for task in self.reference_sha256.keys() & hashes.keys()
            if self.reference_sha256[task] != hashes[task]
        }

    @classmethod
    def from_json(cls, path: Path | str) -> MetricCalibration:
        return cls(**json.loads(Path(path).read_text()))

    def to_json(self, path: Path | str) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2) + "\n")


# Measured by running the three shipped verifiers over `sample_submission/` at the default
# salt on 2026-10-31. That submission scores exactly 1.0 by definition, so its costs *are*
# the reference costs. Reproduce with:
#
#     vqcc verify <sample_submission> --verifiers <verifiers> --write-calibration cal.json
#
# For the record, the estimates these replaced -- back-derived from third-party task ratios
# in Kaggle discussion topic 754174 -- came in 0.02%, 0.006% and 0.04% low respectively, so
# the public frontier's self-reported numbers held up.
DEFAULT_CALIBRATION = MetricCalibration(
    reference_costs={
        "invsqrt": 2_076_774.233715282,
        "dqi": 24_439_463_583_672.93,
        "dfthcsos": 3_355_599_588.3155365,
    },
    provenance=(
        "measured: official verifiers over sample_submission/ at the default salt, "
        "2026-10-31 (200,540,928 / 102,400 / 254,464 verified shots)"
    ),
    measured=True,
    reference_sha256={
        "invsqrt": "e2bf8eb64350a4460a271fc7f1d591ef125b1e276377e9a0f64b96ab3ca321e1",
        "dqi": "959471d05a26373b4f14ca0cae2d8639ea25f18170d9f6532139aa51ccc04995",
        "dfthcsos": "276d464ff14da828b6b33e8487d99975494cbdf0c3cb2bd37dd3c3ddc4002695",
    },
)


def rotation_t_count(eps: float) -> float:
    """T gates to synthesise one arbitrary-angle rotation to operator-norm error `eps`."""
    if not 0 < eps < 1:
        raise ValueError(f"synthesis error must be in (0, 1), got {eps}")
    return ROTATION_T_SLOPE * math.log2(1 / eps) + ROTATION_T_OFFSET


def toffoli_equivalents(
    toffoli_count: float, t_count: float, rotation_count: float, rotation_eps: float | None
) -> float:
    """The invsqrt verifier's `C = Toffolis + T/4`, with rotations priced into T.

    `toffoli_count` and `t_count` are per-shot means from the verifier, not static counts.
    """
    total_t = t_count
    if rotation_count:
        if rotation_eps is None:
            raise ValueError(
                f"{rotation_count} rotations need the synthesis error the verifier chose; "
                "it depends on the accuracy margin the simulated error leaves over"
            )
        total_t += rotation_count * rotation_t_count(rotation_eps)
    return toffoli_count + total_t / T_PER_TOFFOLI_EQUIVALENT


def clopper_pearson_upper_bound(
    failures: int, shots: int, log2_delta: float = DEFAULT_LOG2_DELTA
) -> float:
    """One-sided Clopper-Pearson upper bound at confidence `1 - 2^log2_delta`.

    The same bound the dqi verifier scores on. Zero failures in 102,400 shots still bounds
    the rate near 8.7e-4, which is why a flawless decoder pays R = 1.19.
    """
    if shots <= 0 or not 0 <= failures <= shots:
        raise ValueError(f"need 0 <= failures <= shots and shots > 0, got {failures}/{shots}")
    if not (math.isfinite(log2_delta) and log2_delta <= 0.0):
        raise ValueError(f"log2_delta must be finite and <= 0, got {log2_delta}")
    if failures == shots:
        return 1.0
    if log2_delta == 0.0:
        return failures / shots
    # The exact bound is the Beta inverse survival function. `isf` takes the tail
    # probability directly, so delta = 2^-128 never has to be formed as 1 - delta.
    return float(beta.isf(2.0**log2_delta, failures + 1, shots - failures))


def dqi_repetition_factor(failure_rate: float) -> float:
    """R(eps) = ln 2 / -ln(1/2 + 2 sqrt(eps)); R(0) = 1 and R diverges at eps = 1/16."""
    if not 0 <= failure_rate < MAX_DQI_FAILURE_RATE:
        raise ValueError(
            f"failure rate {failure_rate} outside [0, {MAX_DQI_FAILURE_RATE}); "
            "R is non-positive there and the verifier rejects the submission"
        )
    return math.log(2) / -math.log(0.5 + 2 * math.sqrt(failure_rate))


def invsqrt_cost(num_qubits: int, toffoli_equivalent_count: float) -> float:
    return num_qubits * toffoli_equivalent_count


def dqi_cost(num_qubits: int, toffoli_count: float, eps_plus: float) -> float:
    """`toffoli_count` is the per-shot mean; `eps_plus` is the Clopper-Pearson bound."""
    return dqi_repetition_factor(eps_plus) * num_qubits**2 * toffoli_count


def lambda_eff(e_gap_ha: float, lcu_lambda: float) -> float:
    """lambda_eff = sqrt(E_gap (2 Lambda - E_gap)), the chemistry task's leading factor."""
    inner = e_gap_ha * (2 * lcu_lambda - e_gap_ha)
    if inner < 0:
        raise ValueError(
            f"E_gap={e_gap_ha} and Lambda={lcu_lambda} give a negative "
            "E_gap (2 Lambda - E_gap); lambda_eff is not real there"
        )
    return math.sqrt(inner)


def min_lcu_lambda(c0: float) -> float:
    """The smallest `Lambda` a chemistry submission can reach, over both floors on it.

    `MIN_LAMBDA` is the verifier's rejection rule; this is the stronger, derived bound.
    `D = 2 Lambda - E_gap` rearranges to `2P + 2(lambda_h2 - sum(W^2)/4) - c0`, in which the
    negative eigenvalue mass cancels. Both remaining terms are non-negative -- `P` is a mass,
    and `lambda_h2 - sum(W^2)/4` sums `((|w|_1 + |W + sum w|)^2 - W^2) / 4` terms, each
    non-negative because `|w|_1 + |W + sum w| >= |W|`. So `D >= -c0`, and `Lambda` is smallest
    when `D` sits at that bound and `E_gap` at its floor.

    `c0` is fixed by the reference Hamiltonian and the frozen shift, so no factorization can
    cross this. For FeMoCo the derived bound dominates by 5.7x, which is why the one-norm has
    far less room than `MIN_LAMBDA` alone suggests.
    """
    return max(MIN_LAMBDA, (MIN_E_GAP_HA - c0) / 2.0)


def lambda_eff_floor(c0: float) -> float:
    """`lambda_eff` at `min_lcu_lambda(c0)` with `E_gap` at its floor.

    Unreachable rather than merely a bound: it needs `P = 0` and `lambda_h2 = sum(W^2)/4`, the
    second forcing `w = 0`, which destroys the two-body residual. It bounds the search anyway.
    """
    return lambda_eff(MIN_E_GAP_HA, min_lcu_lambda(c0))


def dfthcsos_cost(
    num_qubits: int, toffolis: float, e_gap_ha: float, lcu_lambda: float
) -> float:
    """Q and T cover the whole block encoding, not only the submitted ANGLE circuits.

    The verifier adds the ANGLE family's simulated Toffoli mean to a `heuristic` part it
    derives analytically from the problem shape (orbitals, R, B, C, `brot`), and takes Q as
    the shared persistent-plus-temporary workspace over both.
    """
    return lambda_eff(e_gap_ha, lcu_lambda) * toffolis * num_qubits


def cost_from_report(task: str, report: Mapping[str, Any]) -> float:
    """That task's cost, read from the fields its verifier reports.

    The formulas and the names of their inputs belong together: this is what the model says
    the cost is, given a verifier's own output. Comparing it against that output's `score`
    is how the transcription stays honest -- see `verify.Verification.model_cost`.
    """
    match require_task(task):
        case "invsqrt":
            return invsqrt_cost(
                report["num_qubits"],
                toffoli_equivalents(
                    report["toffoli_count"],
                    report["t_count"],
                    report["rotation_count"],
                    report.get("rotation_eps"),
                ),
            )
        case "dqi":
            return dqi_cost(
                report["num_qubits"],
                report["toffoli_count"],
                report["effective_error_rate"],
            )
        case _:
            # Rebuild lambda_eff from E_gap and Lambda rather than reusing the reported
            # value, so the identity behind the chemistry cost is checked too.
            metrics = report["metrics"]
            return dfthcsos_cost(
                report["num_qubits"],
                report["toffoli_count"],
                metrics["E_gap"],
                metrics["Lambda"],
            )


# The two dicts that carry these terms name the qubit ones differently, and that difference is
# a trap worth stating: in `heuristic_costs`' output `persistent_qubits` and `temporary_qubits`
# are the *heuristic* figures, while in a verifier report the same two keys are the totals over
# both families and the heuristic ones carry a `heuristic_` prefix. Reading a report through the
# closed-form names silently picks up the totals, so each source gets its own mapping.
HEURISTIC_FIELDS = {
    "outer_toffoli": "outer_state_prep_toffoli",
    "inner_toffoli": "inner_state_prep_toffoli",
    "givens_toffoli": "givens_rotations_toffoli",
    "heuristic_persistent_qubits": "persistent_qubits",
    "heuristic_temporary_qubits": "temporary_qubits",
}
REPORT_FIELDS = HEURISTIC_FIELDS | {
    "heuristic_persistent_qubits": "heuristic_persistent_qubits",
    "heuristic_temporary_qubits": "heuristic_temporary_qubits",
}


@dataclass(frozen=True)
class ShapeResources:
    """`T_total` and `Q_total`, split into what a submission buys and what it compiles free.

    The split is the point. Every heuristic term follows analytically from the shape
    `(R, B, C)` and `brot`, so a submission moves it only by refitting the factorization --
    which it pays for in accuracy. The ANGLE terms are measurements of circuits the submitter
    writes, so cheaper loaders for the *same* words score better at *identical* accuracy.
    """

    outer_toffoli: int
    inner_toffoli: int
    givens_toffoli: int
    heuristic_persistent_qubits: int
    heuristic_temporary_qubits: int
    angle_toffoli: float
    angle_word_qubits: int
    angle_temporary_qubits: int

    @classmethod
    def from_terms(
        cls,
        terms: Mapping[str, Any],
        fields: Mapping[str, str],
        angle_toffoli: float,
        angle_word_qubits: int,
        angle_temporary_qubits: int,
    ) -> ShapeResources:
        """The heuristic half out of `fields`, the ANGLE half measured and passed in."""
        return cls(
            **{name: int(terms[key]) for name, key in fields.items()},
            angle_toffoli=angle_toffoli,
            angle_word_qubits=angle_word_qubits,
            angle_temporary_qubits=angle_temporary_qubits,
        )

    @classmethod
    def from_heuristic(
        cls,
        heuristic: Mapping[str, Any],
        angle_toffoli: float,
        angle_word_qubits: int,
        angle_temporary_qubits: int = 0,
    ) -> ShapeResources:
        """From `dfthcsos.verify.heuristic_costs`' output for a shape, built or not."""
        return cls.from_terms(
            heuristic, HEURISTIC_FIELDS, angle_toffoli, angle_word_qubits, angle_temporary_qubits
        )

    @classmethod
    def from_parts_record(cls, part: Mapping[str, Any]) -> ShapeResources:
        """The lumped split a parts record supports, with no verifier tree on the path.

        A record carries the three sums -- `T_total`, `T_angles` and `Q_total` -- and `J`, but
        not the heuristic qubit split, which only `heuristic_costs` produces. So the heuristic
        half is lumped: its whole Toffoli count into `givens_toffoli` and its whole qubit count
        into `heuristic_persistent_qubits`, both temporaries zero. Every quantity derived from
        here reads the sums rather than the terms, so `toffoli`, `qubits` and `cost_proxy` are
        exact.

        Zero angle temporaries is exact too, under re-batching as well as at the recorded
        batching, but by an argument rather than by construction: folding the heuristic
        workspace into the persistent count is only faithful while that workspace stays the
        larger side of the `max` in `qubits`, and `ANGLE_WIDTH` caps a family's own temporaries
        at `MAX_ANGLE_QUBITS` minus its selection register. That holds for every shape this is
        used on and is pinned against the verifier-backed split by test, not assumed.
        """
        brot, batch = int(part["brot"]), int(part["J"])
        words = batch * brot
        return cls(
            outer_toffoli=0,
            inner_toffoli=0,
            givens_toffoli=int(float(part["T_total"]) - float(part["T_angles"])),
            heuristic_persistent_qubits=int(part["Q_total"]) - words,
            heuristic_temporary_qubits=0,
            angle_toffoli=float(part["T_angles"]),
            angle_word_qubits=words,
            angle_temporary_qubits=0,
        )

    @classmethod
    def from_report(cls, report: Mapping[str, Any]) -> ShapeResources:
        """Both halves from a verifier report, which publishes every term separately.

        The one path needing neither the FCIDUMP nor a cost model, so it is what pins the
        closed-form split against a submission the verifier actually scored.
        """
        resources = report["resources"]
        return cls.from_terms(
            resources,
            REPORT_FIELDS,
            angle_toffoli=float(resources["angles_toffoli"]),
            angle_word_qubits=int(resources["angles_persistent_qubits"]),
            angle_temporary_qubits=int(resources["angles_temporary_qubits"]),
        )

    @property
    def heuristic_toffoli(self) -> int:
        return self.outer_toffoli + self.inner_toffoli + self.givens_toffoli

    @property
    def toffoli(self) -> float:
        return self.angle_toffoli + self.heuristic_toffoli

    @property
    def qubits(self) -> int:
        """The verifier's `Q_total`: the workspace is shared, so it is a max and not a sum."""
        return (
            self.heuristic_persistent_qubits
            + self.angle_word_qubits
            + max(self.heuristic_temporary_qubits, self.angle_temporary_qubits)
        )

    @property
    def cost_proxy(self) -> float:
        """`T * Q`: everything in the chemistry cost that a submission still moves.

        `lambda_eff` is the third factor and divides out of any ratio of two costs, so this
        is what a projected change in resources is measured on.
        """
        return self.toffoli * self.qubits

    @property
    def angle_share(self) -> float:
        """The fraction of `T_total` that costs no accuracy."""
        return self.angle_toffoli / self.toffoli

    def angle_cut_for(self, cost_cut: float) -> float:
        """How much cheaper the ANGLE family must get to take `cost_cut` off `T * Q` alone.

        Above 1.0 the target is out of reach from the free share even with the loaders at zero
        Toffolis, which is the question worth asking before writing a single circuit.
        """
        return self.toffoli * cost_cut / self.angle_toffoli

    def blend_from(self, blend: float, before: ShapeResources) -> float:
        """Where `before`'s `blend` lands on these resources.

        A change in `T * Q` alone moves the blend by the cube root of its inverse, which is
        `blend_factor` on the two cost proxies -- so a hypothetical ANGLE family can be priced
        with no calibration and no verifier run.
        """
        return blend * blend_factor(self.cost_proxy, before.cost_proxy)

    def with_angle_family(
        self, family: AngleFamily, ands_per_circuit: float, temporary_qubits: int
    ) -> ShapeResources:
        """The same heuristic half with a differently batched ANGLE family priced onto it.

        The heuristic terms depend on the factorization shape and `brot`, not on `J`, so
        re-batching the loaders moves only the three ANGLE fields -- which is what makes `J`
        a lever that costs no accuracy. `temporary_qubits` is asked for rather than taken from
        the family's width headroom: a construction need not saturate the cap, and the
        submitted one does not (126 qubits, 75 temporary, against 77 of headroom).
        """
        return replace(
            self,
            angle_toffoli=family.toffoli(ands_per_circuit),
            angle_word_qubits=family.output_bits,
            angle_temporary_qubits=temporary_qubits,
        )


def word_registers(
    qubits_total: int, heuristic_persistent: int, heuristic_temporary: int, brot: int
) -> int:
    """`J`, the angle-word register count, recovered from a recorded `Q_total`.

    The one resource term a verifier report states only as a product
    (`angles_persistent_qubits = J * brot`) and a parts record not at all. Everything else in
    `Q` is closed-form in the shape, so J is whatever is left -- exact, not fitted, and it
    fails loudly rather than rounding if the recorded total does not decompose.

    Takes the two heuristic qubit counts as numbers rather than a dict: the key that holds each
    one depends on which of the two sources it came from, per `HEURISTIC_FIELDS` above.
    """
    words = qubits_total - heuristic_persistent - heuristic_temporary
    if words <= 0 or words % brot:
        raise ValueError(
            f"Q_total {qubits_total} leaves {words} qubits for the word registers, "
            f"which is not a positive multiple of brot={brot}"
        )
    return words // brot


def and_counting_floor(output_bits: int, input_bits: int, num_addresses: int) -> int:
    """The fewest AND gates any ANGLE circuit of this shape can have, for a generic table.

    An ANGLE_k circuit is not a lookup table but a GF(2) circuit: the verifier charges
    Toffoli-like gates (`CCX`, `CCZ`, `CSWAP`) only, and a conditioned gate only on the shots
    whose bit is set, so XOR is free and so is the measurement-based uncompute, whose fixups
    are Clifford. A straight-line circuit with `input_bits` inputs, `A` AND gates and
    `output_bits` outputs realises at most

        2 ** (2*A*n + A**2 + A + b*(n + A + 1))

    distinct maps -- each AND operand is affine in the inputs and the `i - 1` earlier ANDs, and
    each output is affine in the inputs and all `A` ANDs. The circuit has to reproduce a
    specific `b`-bit word on each of `num_addresses` constrained addresses, one of
    `2 ** (b * num_addresses)` tables, so all but a vanishing fraction of tables need

        A**2 + A*(2n + 1 + b) + b*(n + 1) >= b * num_addresses.

    A counting argument, so it bounds almost every table rather than proving anything about
    one: a table with structure can fall far below it.
    """
    if min(output_bits, input_bits, num_addresses) <= 0:
        raise ValueError(
            f"need positive bit and address counts, got {output_bits}, "
            f"{input_bits}, {num_addresses}"
        )
    linear = 2 * input_bits + 1 + output_bits
    free = output_bits * (num_addresses - input_bits - 1)
    return math.ceil((-linear + math.sqrt(linear**2 + 4 * max(free, 0))) / 2)


@dataclass(frozen=True)
class AngleFamily:
    """One admissible ANGLE batching: `K = ceil(75 / J)` circuits of `J * brot` output bits.

    `J` is a free parameter of a submission -- the contract fixes only that the family covers
    all `NUM_ANGLE_WORDS` angles and that no circuit exceeds `MAX_ANGLE_QUBITS`. Raising it
    trades circuits against width: `K` falls as `1 / J` while a circuit's AND floor grows only
    as `sqrt(J)`, so the family gets cheaper in Toffolis and dearer in qubits.
    """

    batch_size: int
    brot: int
    num_addresses: int

    def __post_init__(self) -> None:
        if min(self.batch_size, self.brot, self.num_addresses) <= 0:
            raise ValueError(f"batch size, brot and address count must be positive, got {self}")

    @classmethod
    def from_resources(
        cls, resources: ShapeResources, brot: int, num_addresses: int
    ) -> AngleFamily:
        """The batching a measured submission used, with `J` read off its word registers."""
        words = resources.angle_word_qubits
        if words <= 0 or words % brot:
            raise ValueError(
                f"{words} angle-word qubits is not a positive multiple of brot={brot}"
            )
        return cls(words // brot, brot, num_addresses)

    @property
    def address_bits(self) -> int:
        """The selection width the verifier demands: extra bits would be uncharged workspace."""
        return (self.num_addresses - 1).bit_length()

    @property
    def num_circuits(self) -> int:
        return math.ceil(NUM_ANGLE_WORDS / self.batch_size)

    @property
    def output_bits(self) -> int:
        return self.batch_size * self.brot

    @property
    def public_qubits(self) -> int:
        """Selection plus the word registers, the part of a circuit the contract names."""
        return self.address_bits + self.output_bits

    @property
    def free_workspace(self) -> int:
        """Qubits the width cap leaves for AND ancillas."""
        return MAX_ANGLE_QUBITS - self.public_qubits

    @property
    def feasible(self) -> bool:
        """Whether a circuit this wide is admissible at all, by the one width predicate."""
        return not ANGLE_WIDTH.violated(self.public_qubits)

    @property
    def and_floor(self) -> int:
        """The counting floor, granting the circuit its word registers as inputs.

        They arrive holding the previous batch, which is itself a function of the address, so
        they are `output_bits` nonlinear functions available for free -- and the submitted
        circuits do use them, 8 of circuit 13's 516 AND operands being word qubits. A bound
        that withheld them would not be a floor on what a construction may do.
        """
        return and_counting_floor(self.output_bits, self.public_qubits, self.num_addresses)

    @property
    def address_only_floor(self) -> int:
        """The same bound with the word registers withheld -- the weaker claim, kept to name.

        Granting them is what makes `and_floor` a floor on what a construction may do, and
        the difference is large enough to change a verdict: at the submitted shape this reads
        191 against 157, and a sweep against 191 rules out a target the circuits can in fact
        be built toward. Exposed so the page and the tests cite one derivation of it.
        """
        return and_counting_floor(self.output_bits, self.address_bits, self.num_addresses)

    def toffoli(self, ands_per_circuit: float) -> float:
        """The family's charged Toffoli mean: every AND runs once per pass, and there are two."""
        return 2 * self.num_circuits * ands_per_circuit

    def ands_for_toffoli(self, toffoli: float) -> float:
        """The inverse: ANDs per circuit that a family-wide Toffoli mean implies."""
        return toffoli / (2 * self.num_circuits)


def rebatched(
    resources: ShapeResources, family: AngleFamily, ands_per_circuit: float
) -> ShapeResources:
    """`resources` with `family` priced onto it at `ands_per_circuit` AND gates per circuit.

    The temporary count is zero rather than the family's width headroom: `ANGLE_WIDTH` caps a
    circuit at `MAX_ANGLE_QUBITS`, so its workspace cannot exceed that minus its selection
    register, which is below the heuristic workspace of any shape worth sweeping -- it never
    reaches the `max` in `ShapeResources.qubits`. Passing zero states that, where passing the
    headroom would quietly depend on it.
    """
    return resources.with_angle_family(family, ands_per_circuit, 0)


@dataclass(frozen=True)
class AngleFloorRow:
    """One admissible batching, priced at the counting floor and at a measured AND rate."""

    family: AngleFamily
    at_floor: float
    at_rate: float
    ands_for_target: float

    @property
    def batch_size(self) -> int:
        return self.family.batch_size

    @property
    def target_over_floor(self) -> float:
        """What a family reaching the target would cost as a multiple of the bound."""
        return self.ands_for_target / self.family.and_floor

    @property
    def reaches_target(self) -> bool:
        """Whether the target is reachable here at all, loaders at zero Toffolis included."""
        return self.ands_for_target > 0


@dataclass(frozen=True)
class AngleFloorSweep:
    """Every admissible ANGLE batching of one submission, priced without a verifier run.

    `J` is the lever the contract leaves free -- it fixes only that a family covers all
    `NUM_ANGLE_WORDS` angles and that no circuit exceeds `MAX_ANGLE_QUBITS` -- and raising it
    trades circuits against width: `K` falls as `1 / J` while a circuit's AND floor grows only
    as `sqrt(J)`. The width cap ends the sweep.

    Built once and rendered by both `vqcc angle-floor` and the artifact page, so the table a
    reader runs and the chart a reader is shown cannot disagree.
    """

    submitted: AngleFamily
    measured_ands: float
    blend: float
    target: float
    rows: tuple[AngleFloorRow, ...]

    @classmethod
    def over(
        cls,
        resources: ShapeResources,
        brot: int,
        num_addresses: int,
        blend: float,
        target: float = ESCAPE_SCORE,
        max_batch: int = NUM_ANGLE_WORDS,
    ) -> AngleFloorSweep:
        """Sweep `J` on a measured submission, holding its accuracy and its heuristic half."""
        submitted = AngleFamily.from_resources(resources, brot, num_addresses)
        measured = submitted.ands_for_toffoli(resources.angle_toffoli)
        rate = measured / submitted.and_floor
        # The cost proxy that scores `target`; a candidate's own Q then fixes the T it needs.
        wanted = resources.cost_proxy * (blend / target) ** len(TASKS)
        rows = []
        for size in range(1, max_batch + 1):
            family = AngleFamily(size, brot, num_addresses)
            if not family.feasible:
                continue
            floor = family.and_floor
            rows.append(
                AngleFloorRow(
                    family=family,
                    at_floor=rebatched(resources, family, floor).blend_from(blend, resources),
                    at_rate=rebatched(resources, family, rate * floor).blend_from(
                        blend, resources
                    ),
                    ands_for_target=family.ands_for_toffoli(
                        wanted / rebatched(resources, family, 0).qubits
                        - resources.heuristic_toffoli
                    ),
                )
            )
        return cls(submitted, measured, blend, target, tuple(rows))

    @property
    def rate(self) -> float:
        """The submitted loaders' AND count as a multiple of their own counting floor."""
        return self.measured_ands / self.submitted.and_floor

    @property
    def best_floor(self) -> AngleFloorRow:
        return max(self.rows, key=lambda row: row.at_floor)

    @property
    def best_rate(self) -> AngleFloorRow:
        return max(self.rows, key=lambda row: row.at_rate)

    @property
    def permits_target(self) -> bool:
        """Whether any admissible batching could reach the target at the counting floor."""
        return self.best_floor.at_floor >= self.target

    def at(self, batch_size: int) -> AngleFloorRow | None:
        """The row for one batching, or None where the width cap forbids it."""
        return next((row for row in self.rows if row.batch_size == batch_size), None)


def task_ratio(task: str, cost: float, calibration: MetricCalibration) -> float:
    if cost <= 0:
        raise ValueError(f"{task} cost must be positive, got {cost}")
    return calibration.reference_costs[require_task(task)] / cost


def blend_factor(cost: float, reference_cost: float) -> float:
    """What moving one task's cost multiplies the whole blend by.

    The blend is a geometric mean of three ratios, so changing one task's cost scales it by
    the cube root of the inverse change -- whatever the other two costs happen to be. That is
    why 1% off *any* task is worth the same, and it makes a sensitivity computable without a
    calibration or the other two measurements. `score` is the product of these, which is the
    identity rather than a coincidence.
    """
    if cost <= 0 or reference_cost <= 0:
        raise ValueError(f"costs must be positive, got {cost} against {reference_cost}")
    return (reference_cost / cost) ** (1 / len(TASKS))


def score(costs: Mapping[str, float], calibration: MetricCalibration) -> float:
    """Geometric mean of the three per-task improvement ratios."""
    if set(costs) != set(TASKS):
        raise ValueError(f"need a cost for each of {TASKS}, got {sorted(costs)}")
    return math.prod(
        blend_factor(costs[t], calibration.reference_costs[require_task(t)]) for t in TASKS
    )


def cost_cut_for_blend(blend: float, target: float) -> float:
    """The fraction off one task's cost that lifts `blend` to `target`.

    The inverse of `blend_factor`: what a stated score goal costs, in one number, before any
    of it is attributed to lambda_eff, T or Q.
    """
    if blend <= 0 or target <= 0:
        raise ValueError(f"blends must be positive, got {blend} and {target}")
    return 1.0 - (blend / target) ** len(TASKS)


def marginal_blend_per_unit(blend: float, total: float) -> float:
    """What one unit off a resource of size `total` is worth, as blend points.

    The derivative of `blend_factor` at the current point, so the Toffoli price and the qubit
    price come from one place: `Q` being ~38x smaller than `T` is the whole reason a qubit is
    worth ~38 Toffolis.
    """
    if blend <= 0 or total <= 0:
        raise ValueError(f"blend and total must be positive, got {blend} and {total}")
    return blend / (len(TASKS) * total)


# Every hard constraint the dfthcsos verifier rejects on, in one table. `sense` is the whole
# difference between a cap and a floor: the verifier rejects an error or a width *above* its
# cap but rejects E_gap or Lambda at or below its floor, so a cap admits a value sitting
# exactly on it and the two floors do not. `paths` locate the quantity in a verifier report;
# two of them mean the binding one, since the 20 mHa cap applies to each accuracy term
# separately. One `violated` predicate now settles which side of a limit is legal, for both
# the pre-flight width gate and a measured submission's headroom -- the two previously drifted
# on which accuracy term binds, and on whether a value sitting exactly on a limit is legal.
AT_OR_BELOW = "at_or_below"
STRICTLY_ABOVE = "strictly_above"


@dataclass(frozen=True)
class Constraint:
    """One hard limit, and how to read the quantity it bounds off a verifier report."""

    name: str
    limit: float
    sense: Literal["at_or_below", "strictly_above"]
    paths: tuple[tuple[str, str], ...]

    @property
    def legal_at_limit(self) -> bool:
        return self.sense == AT_OR_BELOW

    def value(self, report: Mapping[str, Any]) -> float:
        # No float() cast: the width is an int in the verifier's report and is serialised as
        # one, and `max` over a single path returns that value unchanged.
        return max(report[section][key] for section, key in self.paths)

    def violated(self, value: float) -> bool:
        return value > self.limit if self.legal_at_limit else value <= self.limit

    def fraction_of_limit_used(self, value: float) -> float:
        """Closeness to rejection, normalised over both directions: 1.0 sits on the limit."""
        return value / self.limit if self.legal_at_limit else self.limit / value

    def margin(self, report: Mapping[str, Any]) -> dict[str, Any]:
        value = self.value(report)
        return {
            "value": value,
            "limit": self.limit,
            "sense": self.sense,
            "fraction_of_limit_used": self.fraction_of_limit_used(value),
        }


ANGLE_WIDTH = Constraint(
    "angle_qubits", MAX_ANGLE_QUBITS, AT_OR_BELOW, (("resources", "angles_qubits"),)
)
DFTHCSOS_CONSTRAINTS = (
    Constraint(
        "hamiltonian_error_ha",
        MAX_HAMILTONIAN_ERROR_HA,
        AT_OR_BELOW,
        (("metrics", "error_1"), ("metrics", "error_2")),
    ),
    ANGLE_WIDTH,
    Constraint("e_gap_ha", MIN_E_GAP_HA, STRICTLY_ABOVE, (("metrics", "E_gap"),)),
    Constraint("lcu_lambda", MIN_LAMBDA, STRICTLY_ABOVE, (("metrics", "Lambda"),)),
)


def dfthcsos_rejections(angle_circuit_qubits: Sequence[int]) -> list[str]:
    """Why the chemistry submission would be rejected, from the ANGLE circuit widths alone.

    The width is the only one of the four constraints knowable before simulating: the two
    accuracy terms, E_gap and Lambda all come out of a verifier run, and
    `dfthcsos_constraint_margins` checks those against the same table once a report exists.
    The value here is catching a doomed candidate from the widths, before an hour of it.
    """
    too_wide = [q for q in angle_circuit_qubits if ANGLE_WIDTH.violated(q)]
    if not too_wide:
        return []
    return [f"{len(too_wide)} ANGLE circuit(s) wider than {MAX_ANGLE_QUBITS} qubits: {too_wide}"]


def lcu_vectors(rank: int, bases: int) -> int:
    """How many LCU vectors a factorization of this shape makes the ANGLE circuits encode.

    `verify.py` builds `num_orbitals + R*B` of them, which is the quantity the angle circuits
    carry and the one a cheaper shape has to reduce. It is checkable against a submission
    without simulating: it is the row count of the npz's `angle_words`.
    """
    return NUM_ORBITALS + rank * bases


def dfthcsos_constraint_margins(report: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """How close a measured chemistry submission sits to each constraint that would reject it."""
    return {c.name: c.margin(report) for c in DFTHCSOS_CONSTRAINTS}
