"""Generators behind one interface (`core.Generator`).

`MockGenerator` runs anywhere and lets the verifier be built/tested without a GPU.
`RFdiffusionAdapter` / `BoltzGenAdapter` are the real generators — they run on the
A100 and parse their output PDB into `CoordinationSite`s; not runnable on this
machine, so they raise until wired during the A100 setup step (see docs/specs).
"""

from __future__ import annotations

import glob
from dataclasses import replace
from pathlib import Path

import numpy as np

from .core import BinderDesign, CoordinationSite
from .geometry.parse import coordination_site_from_pdb


def load_designs(
    glob_pattern: str,
    pdb_element: str = "NI",
    metal_label: str = "Ni2+",
    generator: str = "design",
    cutoff: float = 2.8,
) -> list[BinderDesign]:
    """Parse design PDBs matching a glob into BinderDesigns, recording each `source`
    path and skipping any file with no metal site. The shared ingestion any generator
    adapter or analysis script uses."""
    designs = []
    for pdb in sorted(glob.glob(glob_pattern, recursive=True)):
        try:
            site = coordination_site_from_pdb(pdb, pdb_element, metal_label, cutoff)
        except ValueError:
            continue
        designs.append(BinderDesign(Path(pdb).stem, site, generator, float("nan"), source=pdb))
    return designs

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
            jittered = replace(base, ligand_xyz=base.ligand_xyz + rng.normal(0, 0.05, size=base.ligand_xyz.shape))
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
    """Primary POC generator — ingests RFdiffusionAA output PDBs into BinderDesigns.

    RFdiffusionAA is launched out-of-band on a GPU box via apptainer; this adapter
    parses the resulting design PDBs into CoordinationSites. Like every adapter it is
    generator-blind on the far side — the verifier only ever sees the site.
    """

    def __init__(
        self,
        output_dir: str | Path,
        pdb_element: str = "NI",
        metal_label: str = "Ni2+",
        cutoff: float = 2.8,
    ):
        self.output_dir = Path(output_dir)
        self.pdb_element = pdb_element
        self.metal_label = metal_label
        self.cutoff = cutoff

    def design(self, target: str | None = None, n: int | None = None) -> list[BinderDesign]:
        designs = load_designs(str(self.output_dir / "*.pdb"), self.pdb_element,
                               target or self.metal_label, "rfdiffusion_aa", self.cutoff)
        if not designs:
            raise FileNotFoundError(f"no design PDBs in {self.output_dir}")
        return designs[:n]


class BoltzGenAdapter:
    """Second real generator — proves the verifier is generator-blind."""

    def design(self, target: str, n: int = 5) -> list[BinderDesign]:
        raise NotImplementedError("BoltzGen wired in a later A100 step")
