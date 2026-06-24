"""Generators behind one interface (`core.Generator`).

`MockGenerator` runs anywhere and lets the verifier be built/tested without a GPU.
`RFdiffusionAdapter` / `BoltzGenAdapter` are the real generators — they run on the
A100 and parse their output PDB into `CoordinationSite`s; not runnable on this
machine, so they raise until wired during the A100 setup step (see docs/specs).
"""

from __future__ import annotations

import numpy as np

from .core import BinderDesign, CoordinationSite

# Idealised octahedral ligand directions (unit vectors along ±x, ±y, ±z).
_OCTAHEDRON = np.array(
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float
)


def octahedral_site(metal: str, bond: float = 2.09, center=(0.0, 0.0, 0.0)) -> CoordinationSite:
    c = np.asarray(center, dtype=float)
    return CoordinationSite(
        metal=metal,
        metal_xyz=c,
        ligand_xyz=c + _OCTAHEDRON * bond,
        ligand_elems=("N", "N", "O", "O", "N", "O"),
    )


class MockGenerator:
    """Deterministic stand-in generator: near-ideal octahedral sites with small jitter."""

    def __init__(self, seed: int = 0):
        self.seed = seed

    def design(self, target: str, n: int = 5) -> list[BinderDesign]:
        rng = np.random.default_rng(self.seed)
        designs = []
        for i in range(n):
            base = octahedral_site(target)
            jittered = CoordinationSite(
                metal=base.metal,
                metal_xyz=base.metal_xyz,
                ligand_xyz=base.ligand_xyz + rng.normal(0, 0.05, size=base.ligand_xyz.shape),
                ligand_elems=base.ligand_elems,
            )
            designs.append(
                BinderDesign(
                    sequence=f"MOCK{i}",
                    site=jittered,
                    generator="mock",
                    generator_confidence=float(rng.uniform(0.5, 0.9)),
                )
            )
        return designs


class RFdiffusionAdapter:
    """Primary POC generator: RFdiffusionAA on the A100, output PDB → CoordinationSite."""

    def design(self, target: str, n: int = 5) -> list[BinderDesign]:
        raise NotImplementedError("RFdiffusionAA runs on pi-a100-80gb; wired during A100 setup")


class BoltzGenAdapter:
    """Second real generator + on-site sponsor tool — proves verifier generator-blindness."""

    def design(self, target: str, n: int = 5) -> list[BinderDesign]:
        raise NotImplementedError("BoltzGen wired at the event / later A100 step")
