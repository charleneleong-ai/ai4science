"""Reference coordination geometry — the pluggable oracle.

`MockReference` carries hand-set values (a stand-in for the pre-event PDB pull).
`PDBReference` / a CSD-Mogul provider drop in behind the same
`ReferenceDistribution` interface — the verifier never changes.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

_PDB_DATA = Path(__file__).parent.parent / "data" / "pdb_reference.json"


@dataclass(frozen=True)
class MetalGeometry:
    """Empirical coordination geometry for a metal ion.

    `coordination_number` is the modal CN (drives the score); `cn_range` is the
    observed (min, max) CN — a site within it is trusted, so the reference (which
    measured the spread), not the verifier, owns the tolerance.
    """

    metal: str
    coordination_number: int
    bond_length_mean: float  # Angstrom
    bond_length_std: float
    cn_range: tuple[int, int] = (1, 12)


class ReferenceDistribution(Protocol):
    def geometry(self, metal: str) -> MetalGeometry: ...


class MockReference:
    """Hand-set reference values — the stand-in for the PDB/CSD pull.

    Ni2+: octahedral, CN6, Ni–(N/O) ≈ 2.09 Å; Cu2+: ~5-coordinate, ≈ 2.00 Å.
    Swapped for `PDBReference`, then CSD/Mogul, behind the same interface.
    """

    _TABLE = {
        "Ni2+": MetalGeometry("Ni2+", 6, 2.09, 0.08, cn_range=(5, 7)),
        "Cu2+": MetalGeometry("Cu2+", 5, 2.00, 0.12, cn_range=(4, 6)),
    }

    def geometry(self, metal: str) -> MetalGeometry:
        if metal not in self._TABLE:
            raise KeyError(f"no reference geometry for {metal!r}")
        return self._TABLE[metal]


class PDBReference:
    """Empirical geometry from a pull of metal sites from the public PDB.

    Loads the cached distribution built by scripts/build_pdb_reference.py (mean/std
    of real first-shell metal–donor bonds + modal coordination number). Same
    interface as MockReference, so the verifier is unchanged. CSD/Mogul drops in
    here next, behind the same `geometry()`.
    """

    def __init__(self, path: str | Path = _PDB_DATA):
        self._table = {
            metal: MetalGeometry(
                g["metal"], g["coordination_number"], g["bond_length_mean"],
                g["bond_length_std"], cn_range=tuple(g["cn_range"]),
            )
            for metal, g in json.loads(Path(path).read_text()).items()
        }

    def geometry(self, metal: str) -> MetalGeometry:
        if metal not in self._table:
            raise KeyError(f"no PDB reference geometry for {metal!r}")
        return self._table[metal]
