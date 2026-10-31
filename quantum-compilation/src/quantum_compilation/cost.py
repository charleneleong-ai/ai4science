"""The competition's three cost functions, and the calibration they depend on.

Transcribed from the Evaluation page (https://www.kaggle.com/competitions/
google-verifiable-quantum-circuit-compilation/overview/evaluation) and the published metric
at https://www.kaggle.com/code/metric/google-quantum-ai-metric:

    Score = (prod_b  cost_ref_b / cost_b) ** (1/3)

    invsqrt   cost = Q * C                with C = Toffolis + T/4
    dqi       cost = R(eps_plus) * Q^2 * T
    dfthcsos  cost = lambda_eff * T * Q   with lambda_eff = sqrt(E_gap (2 Lambda - E_gap))

Two numbers the metric owns and has not published sit in `MetricCalibration`: each task's
reference cost, and the price in T gates of a rotation Clifford+T cannot express exactly.
Both are the same kind of unknown, so they travel together in one loadable object with its
own provenance -- and both arrive by filling that object in, not by editing this module.
The verifiers shipped with the competition data remain the authoritative specification.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from quantum_compilation.resources import Resources

TASKS = ("invsqrt", "dqi", "dfthcsos")

# C = Toffolis + T/4, from the invsqrt cost definition.
T_PER_TOFFOLI_EQUIVALENT = 4

# Hard constraints from the data description. A submission is rejected outright, not merely
# scored badly, when one of these is violated. All four are dfthcsos-only.
MAX_ANGLE_QUBITS = 128
MAX_HAMILTONIAN_ERROR_HA = 20e-3
MIN_LAMBDA = 49.0
MIN_E_GAP_HA = 3.914898

# R(eps) is positive only for eps < 1/16: at eps = 1/16 the denominator -ln(1/2 + 2 sqrt(eps))
# hits zero, so a decoder failing more than 6.25% of shots has unbounded cost.
MAX_DQI_FAILURE_RATE = 1 / 16


class UnpricedRotationError(ValueError):
    """A circuit holds rotations whose T price the calibration does not carry."""


@dataclass(frozen=True)
class MetricCalibration:
    """The metric's own constants, with where they came from.

    `rotation_t_cost` maps a Clifford hierarchy level to its price in T gates. An empty
    mapping means the price is still unknown, and pricing a circuit that needs it raises
    rather than guessing -- a wrong default would silently mis-rank every candidate.
    """

    reference_costs: dict[str, float]
    rotation_t_cost: dict[int, float]
    provenance: str
    measured: bool

    def __post_init__(self) -> None:
        if set(self.reference_costs) != set(TASKS):
            raise ValueError(
                f"reference costs must cover exactly {TASKS}, "
                f"got {sorted(self.reference_costs)}"
            )
        bad = {t: c for t, c in self.reference_costs.items() if c <= 0}
        if bad:
            raise ValueError(f"reference costs must be positive, got {bad}")

    @classmethod
    def from_json(cls, path: Path | str) -> MetricCalibration:
        raw = json.loads(Path(path).read_text())
        return cls(
            reference_costs=raw["reference_costs"],
            # JSON object keys are strings; Clifford hierarchy levels are ints.
            rotation_t_cost={int(k): v for k, v in raw["rotation_t_cost"].items()},
            provenance=raw["provenance"],
            measured=raw["measured"],
        )

    def to_json(self, path: Path | str) -> None:
        Path(path).write_text(
            json.dumps(
                {
                    "reference_costs": self.reference_costs,
                    "rotation_t_cost": {
                        str(k): v for k, v in self.rotation_t_cost.items()
                    },
                    "provenance": self.provenance,
                    "measured": self.measured,
                },
                indent=2,
            )
            + "\n"
        )


# Back-derived by multiplying the public frontier's cost by its reported ratio in
# manwithacat's 2026-10-30 briefing (topic 754174) -- third-party, re-verified by that
# author rather than by us, and quoted to four significant figures, so these carry about
# 0.1% of slop. Replace by measuring the competition's own `sample_submission/`, which
# scores exactly 1.0: `vqcc calibrate <sample_submission>` writes the measured file.
PROVISIONAL_CALIBRATION = MetricCalibration(
    reference_costs={
        "invsqrt": 20.57 * 100_941,
        "dqi": 60.59 * 4.0334e11,
        "dfthcsos": 9.93 * 3.3780e8,
    },
    rotation_t_cost={},
    provenance=(
        "provisional: back-derived from third-party task ratios reported in Kaggle "
        "discussion topic 754174 on 2026-10-30; rotation T price unknown"
    ),
    measured=False,
)


def t_count(resources: Resources, calibration: MetricCalibration) -> float:
    """T gates, counting a T-level `Z_POW` as one and pricing finer rotations."""
    missing = sorted(
        level
        for level, _ in resources.fine_rotations
        if level not in calibration.rotation_t_cost
    )
    if missing:
        raise UnpricedRotationError(
            f"{resources.num_fine_rotations} rotation(s) at Clifford hierarchy level(s) "
            f"{missing} have no T price in this calibration ({calibration.provenance})"
        )
    return resources.num_t_like + sum(
        n * calibration.rotation_t_cost[level] for level, n in resources.fine_rotations
    )


def toffoli_equivalents(resources: Resources, calibration: MetricCalibration) -> float:
    """The metric's `C = Toffolis + T gates / 4`."""
    return (
        resources.num_toffoli_like
        + t_count(resources, calibration) / T_PER_TOFFOLI_EQUIVALENT
    )


def invsqrt_cost(num_qubits: int, toffoli_equivalent_count: float) -> float:
    return num_qubits * toffoli_equivalent_count


def dqi_repetition_factor(failure_rate: float) -> float:
    """R(eps) = ln 2 / -ln(1/2 + 2 sqrt(eps)); R(0) = 1 and R diverges at eps = 1/16."""
    if not 0 <= failure_rate < MAX_DQI_FAILURE_RATE:
        raise ValueError(
            f"failure rate {failure_rate} outside [0, {MAX_DQI_FAILURE_RATE}); "
            "R is non-positive there and the cost is unbounded"
        )
    return math.log(2) / -math.log(0.5 + 2 * math.sqrt(failure_rate))


def dqi_cost(num_qubits: int, toffolis: float, failure_rate: float) -> float:
    return dqi_repetition_factor(failure_rate) * num_qubits**2 * toffolis


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
    """Q and T cover the whole block encoding, not only the submitted ANGLE circuits."""
    return lambda_eff(e_gap_ha, lcu_lambda) * toffolis * num_qubits


def task_ratio(task: str, cost: float, calibration: MetricCalibration) -> float:
    if task not in TASKS:
        raise KeyError(f"unknown task {task!r}; expected one of {TASKS}")
    if cost <= 0:
        raise ValueError(f"{task} cost must be positive, got {cost}")
    return calibration.reference_costs[task] / cost


def score(costs: dict[str, float], calibration: MetricCalibration) -> float:
    """Geometric mean of the three per-task improvement ratios."""
    if set(costs) != set(TASKS):
        raise ValueError(f"need a cost for each of {TASKS}, got {sorted(costs)}")
    return math.prod(task_ratio(t, costs[t], calibration) for t in TASKS) ** (
        1 / len(TASKS)
    )


def dfthcsos_rejections(
    angle_circuit_qubits: list[int],
    hamiltonian_error_ha: float | None = None,
    e_gap_ha: float | None = None,
    lcu_lambda: float | None = None,
) -> list[str]:
    """Every reason the chemistry submission would be rejected, in the data description's order.

    Each constraint but the width is checked only when its quantity is supplied: the widths
    come from the submitted circuits, while the Hamiltonian error, E_gap and Lambda come
    from the official verifier. A gate that can only check part of its inputs should say so
    by checking that part, not by assuming a value for the rest.
    """
    reasons: list[str] = []
    too_wide = [q for q in angle_circuit_qubits if q > MAX_ANGLE_QUBITS]
    if too_wide:
        reasons.append(
            f"{len(too_wide)} ANGLE circuit(s) wider than {MAX_ANGLE_QUBITS} qubits: {too_wide}"
        )
    if (
        hamiltonian_error_ha is not None
        and hamiltonian_error_ha > MAX_HAMILTONIAN_ERROR_HA
    ):
        reasons.append(
            f"encoded Hamiltonian is {hamiltonian_error_ha * 1e3:.3f} mHa from the "
            f"reference, over the {MAX_HAMILTONIAN_ERROR_HA * 1e3:.0f} mHa budget"
        )
    if lcu_lambda is not None and lcu_lambda <= MIN_LAMBDA:
        reasons.append(
            f"LCU normalization Lambda={lcu_lambda} is not above {MIN_LAMBDA}"
        )
    if e_gap_ha is not None and e_gap_ha <= MIN_E_GAP_HA:
        reasons.append(f"E_gap={e_gap_ha} Ha is not above {MIN_E_GAP_HA} Ha")
    return reasons
