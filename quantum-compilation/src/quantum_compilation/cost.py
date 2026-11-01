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
from dataclasses import asdict, dataclass, field
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

# The orbital count of the FeMoCo instance every dfthcsos submission factorizes, and the
# baseline factorization's shape (R, B, C) -- the rank, the basis count and the ensemble rank.
# Both come from the shipped verifier rather than from its JSON report, which is why reading
# the report alone cannot recover them: `scoreinfra/reference.py` declares
# `BASELINE_SHAPE = (65, 60, 60)` and documents the submitted arrays as `U` (R, B, N) and
# `W` (R, B, C). The sample submission's npz agrees, carrying `U` of (65, 60, 76).
NUM_ORBITALS = 76
BASELINE_FACTORIZATION = (65, 60, 60)

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


def task_ratio(task: str, cost: float, calibration: MetricCalibration) -> float:
    if cost <= 0:
        raise ValueError(f"{task} cost must be positive, got {cost}")
    return calibration.reference_costs[require_task(task)] / cost


def score(costs: Mapping[str, float], calibration: MetricCalibration) -> float:
    """Geometric mean of the three per-task improvement ratios."""
    if set(costs) != set(TASKS):
        raise ValueError(f"need a cost for each of {TASKS}, got {sorted(costs)}")
    return math.prod(task_ratio(t, costs[t], calibration) for t in TASKS) ** (1 / len(TASKS))


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
