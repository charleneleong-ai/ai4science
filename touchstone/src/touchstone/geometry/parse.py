"""Extract a CoordinationSite from a designed structure (PDB or mmCIF).

The bridge from a generator's raw output to something the verifier can judge: find
the metal atom, then the first-shell donor atoms (N/O/S) coordinating it. Pure and
dependency-light — runs anywhere, independent of which generator wrote the structure
or whether it emitted PDB or CIF.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np

from ..core import CoordinationSite

# Protein/ligand atoms that typically donate to a metal centre.
DONOR_ELEMENTS = frozenset({"N", "O", "S"})
# Solvent is not a donor for our purposes. The geometry prior is built from protein donors only
# (a designed structure carries no waters), so counting water oxygens here would compare a site
# against a reference measured on a different donor set — the same domain error, mirrored. Designs
# have no solvent, so this only bites on real/co-folded structures that do.
SOLVENT_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O"})

Atom = tuple[str, str, np.ndarray]  # (element symbol, residue name, xyz)


def pdb_element(line: str) -> str:
    """PDB element symbol (cols 77–78), falling back to the atom name."""
    el = line[76:78].strip()
    if not el:  # some writers leave the element column blank
        el = line[12:16].strip().lstrip("0123456789")[:1]
    return el.upper()


def pdb_atoms(text: str) -> Iterator[Atom]:
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            yield pdb_element(line), line[17:20].strip().upper(), np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )


def cif_atoms(text: str) -> Iterator[Atom]:
    """Minimal mmCIF `_atom_site` loop reader (element + residue + coordinates)."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        cols, j = [], i + 1
        while j < len(lines) and lines[j].strip().startswith("_"):
            cols.append(lines[j].strip())
            j += 1
        names = [c.split(".", 1)[1] if "." in c else c for c in cols]
        if "_atom_site.type_symbol" not in cols:
            i = j
            continue
        ei, xi, yi, zi = (names.index(k) for k in ("type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"))
        ri = names.index("label_comp_id") if "label_comp_id" in names else None  # absent in minimal CIFs
        while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(("#", "_", "loop_")):
            f = lines[j].split()
            if len(f) > max(ei, xi, yi, zi):
                res = f[ri].upper() if ri is not None and len(f) > ri else ""
                yield f[ei].upper(), res, np.array([float(f[xi]), float(f[yi]), float(f[zi])])
            j += 1
        i = j


def site_from_atoms(
    atoms: Iterator[Atom], metal_element: str, metal_label: str, cutoff: float, min_dist: float
) -> CoordinationSite:
    metal_xyz: np.ndarray | None = None
    donors: list[tuple[str, np.ndarray]] = []
    for el, res, xyz in atoms:
        if el == metal_element.upper():
            metal_xyz = xyz
        elif el in DONOR_ELEMENTS and res not in SOLVENT_RESIDUES:
            donors.append((el, xyz))
    if metal_xyz is None:
        raise ValueError(f"no {metal_element!r} atom found")

    shell = [(el, xyz) for el, xyz in donors if min_dist < np.linalg.norm(xyz - metal_xyz) <= cutoff]
    return CoordinationSite(
        metal=metal_label,
        metal_xyz=metal_xyz,
        ligand_xyz=np.array([xyz for _, xyz in shell]) if shell else np.empty((0, 3)),
        ligand_elems=tuple(el for el, _ in shell),
    )


def coordination_site_from_pdb(
    pdb_path: str | Path, pdb_element: str, metal_label: str, cutoff: float = 2.8, min_dist: float = 1.0
) -> CoordinationSite:
    """Build a CoordinationSite for `pdb_element` (e.g. "NI") in a PDB.

    `metal_label` is the verifier's name for the metal (e.g. "Ni2+"); donor atoms in
    the distance window `(min_dist, cutoff]` Angstrom become the first shell. `min_dist`
    is a physical floor — no metal–ligand bond is shorter than ~1 Å — so it also drops
    unplaced placeholder atoms some generators park on top of the metal.
    """
    return site_from_atoms(
        pdb_atoms(Path(pdb_path).read_text()), pdb_element, metal_label, cutoff, min_dist
    )


def coordination_site_from_cif(
    cif_path: str | Path, metal_element: str, metal_label: str, cutoff: float = 2.8, min_dist: float = 1.0
) -> CoordinationSite:
    """As `coordination_site_from_pdb`, for an mmCIF structure (e.g. BoltzGen output)."""
    return site_from_atoms(
        cif_atoms(Path(cif_path).read_text()), metal_element, metal_label, cutoff, min_dist
    )


def coordination_site(
    path: str | Path, metal_element: str, metal_label: str, cutoff: float = 2.8, min_dist: float = 1.0
) -> CoordinationSite:
    """Format-agnostic: dispatch on file extension (.cif → mmCIF, else PDB)."""
    reader = coordination_site_from_cif if str(path).endswith((".cif", ".mmcif")) else coordination_site_from_pdb
    return reader(path, metal_element, metal_label, cutoff, min_dist)
