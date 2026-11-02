"""Shared structure and table helpers for the pMHC scripts.

Five scripts in this directory had written the same gemmi load triple and the same four-line
CSV epilogue, and the copies were already drifting: three of them indexed `rows[0]` and so
raised IndexError on an empty result, after the measurement had been printed to stdout but
before anything was written. The sibling campaigns hold better versions of both
(`tnf_binder.tables.write_rows`, `tfr1/scripts/csvio.py`), but tnf-binder is not installed in
the environment this runs in and tfr1 sits on an unmerged branch, so these live here and the
shared package those three projects want is recorded as a follow-up rather than attempted four
days before a deadline.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, Mapping, Protocol, runtime_checkable

import gemmi
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.SASA import ShrakeRupley

# Hydrogens are absent from these depositions and are excluded anyway, so every geometric
# measurement here is heavy-atom. The repo's other campaigns call a contact 4.5 A heavy-atom;
# this one reports contacts at an explicit cutoff per call instead of carrying a default,
# because the TCR readout is quoted at two cutoffs and a single constant invited confusion.
HYDROGEN = gemmi.Element("H")

BACKBONE_ATOMS = frozenset(("N", "CA", "C", "O", "OXT"))


@runtime_checkable
class HasRow(Protocol):
    @property
    def row(self) -> Mapping[str, object]: ...


def load(path: Path) -> gemmi.Structure:
    """Read a structure, drop solvent and ligands, and set up entities.

    Stripping waters matters beyond tidiness: `reachability.escape_fraction` counts an atom as
    blocking, so a retained water would close an approach channel that solvent in fact occupies.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    return st


def one_letter(residue_name: str) -> str:
    """Uppercase one-letter code, or '' for anything that is not an amino acid.

    Not `gemmi.one_letter_code`, which lowercases non-standard residues rather than dropping
    them: that silently promotes a modified residue such as MSE into a standard letter, and the
    sequences built here become the novelty gate's denominator and the co-fold's chain inputs.
    tfr1 recorded the same trap after a good structure failed as the 'wrong chain'.
    """
    info = gemmi.find_tabulated_residue(residue_name)
    if info is None or not info.is_amino_acid():
        return ""
    return info.one_letter_code.upper()


def chain_sequence(chain: gemmi.Chain) -> str:
    return "".join(one_letter(r.name) for r in chain)


def target_sequences(path: Path) -> frozenset[str]:
    """Every sequence belonging to the target, from prepare_target.py's reference file."""
    loaded = json.loads(path.read_text())
    values = (
        loaded.values()
        if isinstance(loaded, dict)
        else [s for entry in loaded for s in entry["chains"].values()]
    )
    return frozenset(str(v) for v in values)


def binder_chain(
    pose: Path, structure: gemmi.Structure, targets: frozenset[str]
) -> str:
    """The name of the chain whose sequence is none of the target's.

    Identified by sequence, not by chain name or by length. BoltzGen writes the design as chain
    A and shifts the target to B/C/D, while Chai relabels differently again; an earlier version
    of this excluded chain names and so extracted the 10-residue epitope, which then matched our
    own deposition at 100% identity and reported it as a fold hit. tfr1 recorded the same lesson.

    A chain counts as the target's when its sequence is any contiguous part of a reference
    chain, not only when it matches one whole. Poses generated against a truncated receptor
    carry the groove alone -- an early smoke set holds the 182-residue alpha-1/alpha-2 platform
    where the references hold the full 276-residue chain -- and exact matching called that
    platform a design. A design that really is a verbatim slice of the target is then left with
    no candidate chain and raises, which is the right way to fail.

    `pose` is carried only for the error message; the parse belongs to the caller, so that
    callers needing both the chain and its atoms do not read the file twice.
    """
    chains = {ch.name: chain_sequence(ch) for ch in structure[0]}
    candidates = [
        name
        for name, seq in chains.items()
        if seq and not any(seq in reference for reference in targets)
    ]
    if len(candidates) != 1:
        lengths = {k: len(v) for k, v in chains.items()}
        raise ValueError(
            f"{pose}: expected exactly one non-target chain, found {candidates} among {lengths}"
        )
    return candidates[0]


def heavy_coords(chain: gemmi.Chain) -> tuple[np.ndarray, list[int]]:
    """Heavy-atom coordinates of a chain, with each atom's 1-based residue position."""
    coords: list[list[float]] = []
    positions: list[int] = []
    for index, residue in enumerate(chain, start=1):
        for atom in residue:
            if atom.element == HYDROGEN:
                continue
            coords.append(atom.pos.tolist())
            positions.append(index)
    return np.array(coords).reshape(-1, 3), positions


def sasa_model(path: Path, probe_radius: float = 1.40, points: int = 100) -> Model:
    """Parse once and run Shrake-Rupley once, at atom level.

    The callers previously re-parsed and re-computed per question, so one invocation of
    epitope_geometry ran four full SASA passes over a ~4,800-atom complex for four numbers.
    """
    model = MMCIFParser(QUIET=True).get_structure("x", str(path))[0]
    ShrakeRupley(probe_radius=probe_radius, n_points=points).compute(model, level="A")
    return model


def residues_of(model: Model, chain: str) -> list[Residue]:
    """Polymer residues of one chain, hetero flags dropped, in deposition order."""
    return [r for r in model[chain] if r.id[0] == " "]


def chain_area(model: Model, chain: str) -> float:
    return sum(a.sasa for r in residues_of(model, chain) for a in r)


def sidechain_area(residue: Residue, polar_only: bool = False) -> float:
    """Float even when nothing qualifies: a bare sum() returns int 0 for glycine, which then
    formats as `0` rather than `0.0` and churns every committed table."""
    return float(
        sum(
            atom.sasa
            for atom in residue
            if atom.get_id() not in BACKBONE_ATOMS
            and (not polar_only or atom.element in ("N", "O"))
        )
    )


def write_rows(rows: Iterable[HasRow], out: Path, label: str = "rows") -> int:
    """Write dataclass rows to CSV, creating the parent directory, with LF line endings.

    An empty sequence writes nothing and returns 0 rather than raising: these are measurement
    scripts whose findings have already been printed by the time this runs, and a crash here
    loses the run's output to save nothing.
    """
    materialised = list(rows)
    if not materialised:
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(materialised[0].row), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(r.row for r in materialised)
    return len(materialised)
