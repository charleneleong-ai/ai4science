"""Measured resource counts for a kickmix circuit.

Every number here is read off the circuit by kickmix, and kickmix owns the semantics: the
Clifford-hierarchy level of a `Z_POW` and what counts as Toffoli-like or T-like are its
definitions, taken from its own accessors rather than re-derived here.

Nothing in this module knows the competition's cost model, not even the 1/4 that turns T
gates into Toffoli-equivalents -- that lives in `cost.py`, which consumes these counts. Nor
can these counts be turned into a cost: the scored Toffoli figure is a simulated per-shot
mean, so what is measured here bounds a candidate's cost from above and never gives it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import kickmix as km

# A `Z_POW` past the T level is one Clifford+T cannot express exactly, so it needs synthesis.
# The level marks *what needs* synthesising; it does not price it -- the metric charges every
# such rotation at one circuit-wide synthesis error, regardless of its level.
T_LEVEL = 3


@dataclass(frozen=True)
class Resources:
    """Measured width and non-Clifford content of one circuit.

    Deliberately small and hashable: a search loop keeps thousands of these while the
    parsed `km.Circuit` each came from -- up to 1.46 GB for a public dqi decoder -- is
    discarded.
    """

    num_qubits: int
    num_toffoli_like: int
    num_t_like: int
    num_fine_rotations: int

    @classmethod
    def from_circuit(cls, circuit: km.Circuit) -> Resources:
        stats = circuit.operation_statistics()
        by_level = stats.z_pow_clifford_hierarchy_counts
        return cls(
            num_qubits=circuit.num_qubits,
            num_toffoli_like=stats.num_toffoli_like,
            num_t_like=stats.num_t_like,
            num_fine_rotations=sum(n for lvl, n in by_level.items() if lvl > T_LEVEL),
        )

    @classmethod
    def from_file(cls, path: Path | str) -> Resources:
        return cls.from_circuit(km.Circuit.from_file(path))
