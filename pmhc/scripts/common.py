"""Shared structure and table helpers for the pMHC scripts.

Five scripts in this directory had written the same gemmi load triple and the same four-line
CSV epilogue, and the copies were already drifting: three of them indexed `rows[0]` and so
raised IndexError on an empty result, after the measurement had been printed to stdout but
before anything was written. The sibling campaigns hold better versions of both
(`tnf_binder.tables.write_rows`, `tfr1/scripts/csvio.py`), but tnf-binder is not installed in
the environment this runs in and tfr1 sits on an unmerged branch, so these live here and the
shared package those three projects want is recorded as a follow-up rather than attempted four
days before a deadline. Three more belong in that same follow-up: `touchstone.boltz` already reads
Boltz's per-sample confidence files and aggregates them, `tnf_binder.specs.chain_block` is a
byte-level duplicate of this campaign's co-fold spec block, and `tnf_binder.submission` holds the
Proteinbase upload rules -- accepted molecule classes, the 20-row cap, the 10-250 aa range --
which `submission.py` restates because they are portal rules owned by no campaign. That last one
is the most urgent of the four: the portal has already changed its stated format once, glp1r is a
third track against it, and the next change has to be found in two places.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Protocol

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

# The two epitopes, as the challenge states them, and the peptide positions each generation
# variant nominated. Declared here because prepare_target writes the generation specs from them
# and engagement scores the poses against them; when the table lived in both, a hotspot edit in
# one place would have silently rescored the cohort against a different nomination than it was
# generated for.
PEPTIDE_G12D = "VVVGADGVGK"
PEPTIDE_WT = "VVVGAGGVGK"

VARIANTS: dict[str, tuple[int, ...]] = {
    "bulge_core": (5, 6, 7),
    "bulge_wide": (4, 5, 6, 7, 8),
    "mut_proximal": (5, 6),
    "hub": (7, 8),
}


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


def target_references(path: Path) -> dict[str, str]:
    """Named sequences a design must not resemble: the target's chains and the known receptor.

    Reads either shape prepare_target.py writes -- the role-keyed `novelty_references.json`
    (`hla_heavy`, `b2m`, `peptide_g12d`, `peptide_wt`, `tcr_D`, `tcr_E`) or the co-fold
    condition list -- so either file can be passed. The names are the point: a caller wanting
    the epitopes asks for them by role rather than recovering them with a sequence heuristic.
    """
    loaded = json.loads(path.read_text())
    if isinstance(loaded, dict):
        return {str(k): str(v) for k, v in loaded.items()}
    return {
        f"{entry['name']}:{chain}": sequence
        for entry in loaded
        for chain, sequence in entry["chains"].items()
    }


def target_sequences(path: Path) -> frozenset[str]:
    """Every sequence belonging to the target, discarding the roles."""
    return frozenset(target_references(path).values())


def epitopes_from(refs: Mapping[str, str]) -> tuple[str, str]:
    """The G12D and wild-type epitopes, by role rather than by shape.

    An earlier version filtered the reference values for `len(s) == 10 and startswith("VVVGA")`,
    which worked but built its result by iterating a frozenset -- so the order varied per
    interpreter run. Harmless while only membership was tested, and a trap for anyone indexing it.
    """
    try:
        return refs["peptide_g12d"], refs["peptide_wt"]
    except KeyError as missing:
        raise ValueError(f"reference file has no {missing} entry") from missing


def chain_sequences(structure: gemmi.Structure) -> dict[str, str]:
    return {ch.name: chain_sequence(ch) for ch in structure[0]}


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
    chains = chain_sequences(structure)
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


@dataclass(frozen=True)
class ChainRoles:
    """Which chains of a pose are the binder, the peptide, and the rest of the receptor."""

    design: tuple[str, ...]
    peptide: str
    receptor: tuple[str, ...]


def chain_roles(
    structure: gemmi.Structure,
    epitopes: tuple[str, ...],
    *,
    targets: frozenset[str] | None = None,
    design: tuple[str, ...] | None = None,
    pose: Path | None = None,
) -> ChainRoles:
    """Partition a pose's chains into binder / peptide / receptor, in one pass.

    Four separate derivations of "which chain is which" ran over the same pose before this:
    the design by sequence, the peptide by sequence, the receptor as whatever was left, and the
    control's binder hand-typed on the command line. Nothing asserted they agreed.

    `design` names the binder chains when they are known by construction -- the 7OW6 control's
    receptor is two chains, alpha and beta. It cannot be resolved the way a design is: the TCR
    sits in the novelty references on purpose, so that a regurgitated receptor is caught as
    plagiarism, which leaves `binder_chain` with zero candidates on 7OW6. Otherwise pass
    `targets` and the design is the one chain that is none of the target's.
    """
    chains = chain_sequences(structure)
    if design is None:
        if targets is None:
            raise ValueError("chain_roles needs either design chain names or target sequences")
        design = (binder_chain(pose or Path("<pose>"), structure, targets),)
    missing = [name for name in design if name not in chains]
    if missing:
        raise ValueError(f"{pose}: no such chain(s) {missing} among {sorted(chains)}")
    peptide = next((name for name, seq in chains.items() if seq in epitopes), None)
    if peptide is None:
        raise ValueError(f"{pose}: no chain matches either epitope sequence")
    if peptide in design:
        raise ValueError(f"{pose}: chain {peptide} is both the epitope and a binder chain")
    receptor = tuple(n for n in chains if n not in set(design) and n != peptide)
    return ChainRoles(design=tuple(design), peptide=peptide, receptor=receptor)


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


def rows_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    """Every row of a CSV, keyed by one column. A later row with the same key wins."""
    with path.open() as handle:
        return {r[key]: r for r in csv.DictReader(handle)}


def pose_files(poses: Path) -> list[Path]:
    found = sorted(p for p in poses.rglob("*.cif") if p.is_file())
    if not found:
        raise ValueError(f"{poses} holds no .cif poses")
    return found


def variant_of(design: str) -> str:
    """The generation variant a design came from, by name prefix.

    BoltzGen names designs after the spec stem and prepare_target names each spec after its
    variant, so the prefix is reliable for this cohort -- but the coupling is implicit. The
    variant is also authoritative in the pose path (`gen_out/<variant>/...`), which is the
    sounder source and is recorded as a follow-up.
    """
    for name in VARIANTS:
        if design.startswith(name):
            return name
    raise ValueError(f"{design}: no known variant prefix")


def rounded(value: float | None, places: int = 4) -> float | None:
    return None if value is None else round(value, places)


def wrote(out: Path, count: int) -> str:
    """What to echo after write_rows, which writes nothing when there is nothing to write.

    Every caller used to print `wrote <path> (0 rows)` for a file that did not exist.
    """
    return f"wrote {out} ({count} rows)" if count else f"no rows to write; {out} not created"


def write_rows(rows: Iterable[HasRow], out: Path) -> int:
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
