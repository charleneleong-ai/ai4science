"""Extract a CoordinationSite from a designed structure (PDB).

The bridge from a generator's raw output to something the verifier can judge:
find the metal atom, then the first-shell donor atoms (N/O/S) coordinating it.
Pure and dependency-light — runs anywhere, independent of which generator wrote
the PDB.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core import CoordinationSite

# Protein/ligand atoms that typically donate to a metal centre.
DONOR_ELEMENTS = frozenset({"N", "O", "S"})


def _element(line: str) -> str:
    """PDB element symbol (cols 77–78), falling back to the atom name."""
    el = line[76:78].strip()
    if not el:  # some writers leave the element column blank
        el = line[12:16].strip().lstrip("0123456789")[:1]
    return el.upper()


def coordination_site_from_pdb(
    pdb_path: str | Path,
    pdb_element: str,
    metal_label: str,
    cutoff: float = 2.8,
    min_dist: float = 1.0,
) -> CoordinationSite:
    """Build a CoordinationSite for `pdb_element` (e.g. "NI") in a PDB.

    `metal_label` is the verifier's name for the metal (e.g. "Ni2+"); donor atoms
    in the distance window `(min_dist, cutoff]` Angstrom become the first shell.
    `min_dist` is a physical floor — no metal–ligand bond is shorter than ~1 Å — so
    it also drops unplaced placeholder atoms that some generators park on top of the
    metal (e.g. backbone-only output with sidechains at the origin).
    """
    metal_xyz: np.ndarray | None = None
    donors: list[tuple[str, np.ndarray]] = []

    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = _element(line)
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if el == pdb_element.upper():
            metal_xyz = xyz
        elif el in DONOR_ELEMENTS:
            donors.append((el, xyz))

    if metal_xyz is None:
        raise ValueError(f"no {pdb_element!r} atom found in {pdb_path}")

    shell = [
        (el, xyz) for el, xyz in donors if min_dist < np.linalg.norm(xyz - metal_xyz) <= cutoff
    ]
    elems = tuple(el for el, _ in shell)
    ligand_xyz = np.array([xyz for _, xyz in shell]) if shell else np.empty((0, 3))

    return CoordinationSite(
        metal=metal_label, metal_xyz=metal_xyz, ligand_xyz=ligand_xyz, ligand_elems=elems
    )
