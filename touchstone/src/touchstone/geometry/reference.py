"""Reference coordination geometry — the pluggable oracle.

`MockReference` carries hand-set values (a stand-in for the pre-event PDB pull).
`PDBReference` / a CSD-Mogul provider drop in at the event behind the same
`ReferenceDistribution` interface — the verifier never changes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class MetalGeometry:
    """Empirical coordination geometry for a metal ion (mean ± std of bond length)."""

    metal: str
    coordination_number: int
    bond_length_mean: float  # Angstrom
    bond_length_std: float


class ReferenceDistribution(Protocol):
    def geometry(self, metal: str) -> MetalGeometry: ...


class MockReference:
    """Hand-set reference values — the stand-in for the PDB/CSD pull.

    Ni2+: octahedral, CN6, Ni–(N/O) ≈ 2.09 Å; Cu2+: ~5-coordinate, ≈ 2.00 Å.
    Swapped for `PDBReference` (pre-event) then CSD/Mogul (at the event).
    """

    _TABLE = {
        "Ni2+": MetalGeometry("Ni2+", 6, 2.09, 0.08),
        "Cu2+": MetalGeometry("Cu2+", 5, 2.00, 0.12),
    }

    def geometry(self, metal: str) -> MetalGeometry:
        if metal not in self._TABLE:
            raise KeyError(f"no reference geometry for {metal!r}")
        return self._TABLE[metal]


class PDBReference:
    """Empirical geometry from a scripted pull of metal sites from the public PDB.

    Pre-event reference provider: real logic, no license needed. Wired during prep
    (see docs/specs); raises until then so the mock path is the default.
    """

    def geometry(self, metal: str) -> MetalGeometry:
        raise NotImplementedError("PDBReference pull wired during prep; use MockReference for now")
