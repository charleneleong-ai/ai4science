"""Measured resource counts for a kickmix circuit.

Every number here is read off the circuit by kickmix, and kickmix owns the semantics: the
Clifford-hierarchy level of a `Z_POW` and what counts as Toffoli-like or T-like are its
definitions, taken from its own accessors rather than re-derived here.

Nothing in this module knows the competition's cost model. The metric's coefficients -- the
price in T gates of a rotation Clifford+T cannot express exactly, and the 1/4 that turns T
gates into Toffoli-equivalents -- live in `cost.py`, which consumes these counts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import kickmix as km

# Only `Z_POW` past the T level needs Clifford+T synthesis at a cost the metric sets, so
# that is the one level boundary this module has to name.
T_LEVEL = 3


@dataclass(frozen=True)
class Resources:
    """Measured width and non-Clifford content of one circuit.

    Deliberately small and hashable: a search loop keeps thousands of these while the
    parsed `km.Circuit` each came from -- up to 1.46 GB for a public dqi decoder -- is
    discarded. `fine_rotations` is a sorted tuple rather than a dict so `frozen=True` means
    what it says and candidates can be deduped in a set.
    """

    num_qubits: int
    num_toffoli_like: int
    num_t_like: int
    fine_rotations: tuple[tuple[int, int], ...]

    @classmethod
    def from_circuit(cls, circuit: km.Circuit) -> Resources:
        stats = circuit.operation_statistics()
        by_level = stats.z_pow_clifford_hierarchy_counts
        return cls(
            num_qubits=circuit.num_qubits,
            num_toffoli_like=stats.num_toffoli_like,
            num_t_like=stats.num_t_like,
            fine_rotations=tuple(
                sorted((lvl, n) for lvl, n in by_level.items() if lvl > T_LEVEL)
            ),
        )

    @classmethod
    def from_file(cls, path: Path | str) -> Resources:
        return cls.from_circuit(km.Circuit.from_file(path))

    @property
    def num_fine_rotations(self) -> int:
        return sum(n for _, n in self.fine_rotations)
