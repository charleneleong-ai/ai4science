"""Cut the design targets and emit hotspots, with every index checked.

The ectodomain monomer is 639 modelled residues; co-folding that against a binder
for every design is wasteful when the epitope is local, so each candidate epitope
gets its own cut: all residues whose centroid lies within a radius of the patch
centre, kept as contiguous runs so no fragment is a lone floating residue.

Three traps this guards against, all of which have already cost a campaign:

- **A hand-copied conclusion goes stale in silence.** An earlier revision hardcoded
  the three epitope centres with their conservation figures in a comment. Two of
  the three were, in fact, excluded by `target_sites.py`'s own gate, and nothing
  anywhere raised. This version reads the scan's CSV, looks each chosen centre up
  in it, and refuses to cut a target that the scan did not mark usable -- so a
  threshold change that moves the candidate set breaks the build instead of
  shipping stale geometry. The hotspots stay curated (they are a design choice, and
  were weighted towards the receptor carboxylates a histidine was believed to need
  to titrate against -- a belief since measured and refuted, the lever being burial
  on binding; `target_sites.py` records the numbers), but each one is checked to be
  a member of the patch it claims.
- **BoltzGen `binding` indices are 1-based positions in the written chain, not
  UniProt numbers.** The conversion is computed, and then verified by re-reading
  the file that was actually written -- asserting both the residue code and its
  seqid at each converted index. An off-by-N here misdirects every design.
- **A structure that disagrees with UniProt.** Older TfR1 structures carry Ser142
  where UniProt and the assay construct have Gly142; any mismatch raises.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import gemmi
import typer
from target_sites import ACIDIC, centroid, one_letter, read_sequence, truncated_carboxylates

CUT_RADIUS = 22.0
MIN_RUN = 5
RECEPTOR_CHAIN = "R"


@dataclass(frozen=True)
class Epitope:
    centre: int
    hotspots: tuple[int, ...]
    rationale: str


# Chosen from the patch scan; `--patches` is re-read at build time to confirm each
# centre is still usable and each hotspot still belongs to its patch.
EPITOPES: dict[str, Epitope] = {
    "conserved": Epitope(
        centre=529,
        hotspots=(523, 525, 527, 530, 533),
        rationale="most conserved patch clear of both the transferrin site and the glycans",
    ),
    "carboxylate": Epitope(
        centre=605,
        hotspots=(604, 606, 664, 667, 668),
        rationale="richest in carboxylates; chosen on a criterion since refuted, kept as a "
                  "second conserved epitope (83% human/mouse identity)",
    ),
    "apical": Epitope(
        centre=371,
        hotspots=(208, 211, 348, 350, 371),
        rationale="the brief's recommended epitope, carried as the contrast",
    ),
}

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class Patch:
    centre: int
    members: frozenset[int]
    acidic_members: tuple[int, ...]
    conservation: float
    transferrin_distance: float
    glycan_distance: float
    domain: str


def read_patches(path: Path) -> dict[int, Patch]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        int(r["centre"]): Patch(
            centre=int(r["centre"]),
            members=frozenset(int(n) for n in r["residues"].split()),
            acidic_members=tuple(int(n) for n in r["acidic_members"].split()),
            conservation=float(r["conservation"]),
            transferrin_distance=float(r["transferrin_dist"]),
            glycan_distance=float(r["glycan_dist"]),
            domain=r["domain"],
        )
        for r in rows
    }


def amino_acids(chain: gemmi.Chain) -> list[gemmi.Residue]:
    return [r for r in chain if one_letter(r)]


def contiguous_runs(numbers: Collection[int]) -> list[tuple[int, int]]:
    ordered = sorted(numbers)
    if not ordered:
        return []
    runs: list[tuple[int, int]] = []
    start = previous = ordered[0]
    for number in ordered[1:]:
        if number != previous + 1:
            runs.append((start, previous))
            start = number
        previous = number
    runs.append((start, previous))
    return runs


def mismatches_against_uniprot(residues: list[gemmi.Residue], sequence: str) -> list[str]:
    problems = []
    for residue in residues:
        number = residue.seqid.num
        if not 1 <= number <= len(sequence):
            problems.append(f"{residue.name}{number} outside UniProt 1-{len(sequence)}")
        elif one_letter(residue) != sequence[number - 1]:
            problems.append(f"{number}: structure {one_letter(residue)} "
                            f"vs UniProt {sequence[number - 1]}")
    return problems


def cut_runs(centres: dict[int, gemmi.Position], centre_number: int) -> list[tuple[int, int]]:
    centre = centres[centre_number]
    near = [n for n, position in centres.items() if centre.dist(position) <= CUT_RADIUS]
    return [run for run in contiguous_runs(near) if run[1] - run[0] + 1 >= MIN_RUN]


def residues_in(runs: list[tuple[int, int]]) -> list[int]:
    return [n for start, end in runs for n in range(start, end + 1)]


def write_target(residues: list[gemmi.Residue], wanted: Collection[int], path: Path) -> None:
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    chain = gemmi.Chain(RECEPTOR_CHAIN)
    keep = set(wanted)
    for residue in residues:
        if residue.seqid.num in keep:
            chain.add_residue(residue)
    structure[0].add_chain(chain)
    structure.setup_entities()
    structure.write_pdb(str(path))


def positions_in_written_file(path: Path) -> dict[int, int]:
    """UniProt number -> 1-based position, read back off the artefact the generator
    consumes rather than predicted from what we meant to write."""
    written = amino_acids(load_written_chain(path))
    return {residue.seqid.num: index + 1 for index, residue in enumerate(written)}


def load_written_chain(path: Path) -> gemmi.Chain:
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    return structure[0][RECEPTOR_CHAIN]


def read_measured_contacts(path: Path) -> frozenset[int]:
    """Every receptor residue observed in contact with anything, from stage 1.

    A hotspot is legitimate if it is a member of its patch *or* an observed
    contact. For the apical epitope the observed VP1u contacts are the definition
    of the site, and N348 sits just outside a 12 A patch around 371 while being one
    of the nine residues the brief names -- the patch radius is a scanning
    heuristic, the contact is a measurement, and the measurement wins.
    """
    with path.open(newline="") as handle:
        return frozenset(int(row["residue"]) for row in csv.DictReader(handle))


def check_epitope_against_scan(
    name: str, epitope: Epitope, patches: dict[int, Patch], measured: frozenset[int]
) -> Patch:
    patch = patches.get(epitope.centre)
    if patch is None:
        raise typer.BadParameter(
            f"{name}: centre {epitope.centre} is not a usable patch in the scan. Either the "
            f"thresholds moved or the choice is stale -- re-run target_sites.py and re-pick. "
            f"Usable centres: {sorted(patches)}"
        )
    unjustified = [
        h for h in epitope.hotspots if h not in patch.members and h not in measured
    ]
    if unjustified:
        raise typer.BadParameter(
            f"{name}: hotspots {unjustified} are neither members of patch {epitope.centre} "
            f"nor observed contacts in the measured epitopes"
        )
    return patch


def verify_carboxylates(path: Path) -> None:
    """Refuse to ship a target whose acidic sidechains are incomplete.

    Not because a histidine needs a carboxylate to pair with -- that criterion has since been
    measured and refuted -- but because both pKa methods solve for the electrostatic environment
    of the whole complex, and an acidic sidechain missing its charged atoms is scored as
    uncharged. A target truncated this way is a different molecule electrostatically, and the
    linkage measured against it is not the linkage of the real complex. It does not fail loudly:
    it drives every Delta toward zero and quietly reports no coupling. The GLP1R target handed to
    the generator in the previous campaign was truncated at 9 of 48 residues and nothing noticed
    until the switch stage, three stages and a GPU-week later. Checked against the written file,
    not the in-memory copy, for the same reason `verify_positions` is.
    """
    if missing := truncated_carboxylates(amino_acids(load_written_chain(path))):
        raise typer.BadParameter(
            f"{path.name}: {len(missing)} acidic residues are missing their charged atoms "
            f"({', '.join(missing[:6])}). Both pKa methods would score them as uncharged, so "
            f"every linkage measured against this target would be the linkage of a different "
            f"molecule, reported as absence of coupling rather than as a broken target."
        )


def verify_positions(path: Path, hotspots: Collection[int], sequence: str) -> dict[int, int]:
    positions = positions_in_written_file(path)
    written = amino_acids(load_written_chain(path))
    for number in hotspots:
        if number not in positions:
            raise typer.BadParameter(f"{path.name}: hotspot {number} is not in the written chain")
        residue = written[positions[number] - 1]
        if residue.seqid.num != number or one_letter(residue) != sequence[number - 1]:
            raise typer.BadParameter(
                f"{path.name}: position {positions[number]} holds "
                f"{residue.name}{residue.seqid.num}, expected {sequence[number - 1]}{number}"
            )
    return positions


def build_entry(
    name: str,
    epitope: Epitope,
    patch: Patch,
    residues: list[gemmi.Residue],
    centres: dict[int, gemmi.Position],
    sequence: str,
    out_dir: Path,
    measured: frozenset[int],
) -> dict[str, object]:
    runs = cut_runs(centres, epitope.centre)
    kept = residues_in(runs)
    path = out_dir / f"{name}.pdb"
    write_target(residues, kept, path)
    positions = verify_positions(path, epitope.hotspots, sequence)
    verify_carboxylates(path)

    for number in epitope.hotspots:
        basis = "patch" if number in patch.members else "measured contact"
        typer.echo(f"  {name:11s} hotspot {sequence[number - 1]}{number} "
                   f"-> chain position {positions[number]:3d}  ({basis})")
    typer.echo(f"{name:11s} {len(kept):3d} residues in {len(runs)} runs, "
               f"patch conservation {patch.conservation:.0%}, "
               f"{len(patch.acidic_members)} carboxylates -> {path.name}")

    return {
        "centre": epitope.centre,
        "rationale": epitope.rationale,
        "domain": patch.domain,
        "conservation": patch.conservation,
        "transferrin_distance": patch.transferrin_distance,
        "glycan_distance": patch.glycan_distance,
        "patch_carboxylates": list(patch.acidic_members),
        "residues": len(kept),
        "runs": [f"{a}-{b}" for a, b in runs],
        "hotspots_uniprot": list(epitope.hotspots),
        "hotspots_positions": [positions[h] for h in epitope.hotspots],
        "binding_string": ",".join(str(positions[h]) for h in epitope.hotspots),
        "file": path.name,
    }


@app.command()
def main(
    structure_path: Annotated[Path, typer.Option("--structure", help="11QC mmCIF")],
    human_json: Annotated[Path, typer.Option(help="UniProt P02786 JSON")],
    patches: Annotated[Path, typer.Option(help="surface_patches.csv from target_sites.py")],
    contacts: Annotated[Path, typer.Option(help="epitope_contacts.csv from target_sites.py")],
    out_dir: Annotated[Path, typer.Option(help="where targets and hotspots are written")],
) -> None:
    sequence = read_sequence(human_json)
    structure = gemmi.read_structure(str(structure_path))
    structure.setup_entities()
    residues = amino_acids(structure[0]["A"])

    problems = mismatches_against_uniprot(residues, sequence)
    if problems:
        raise typer.BadParameter(
            f"{structure_path.name} disagrees with UniProt at {len(problems)} positions: "
            + "; ".join(problems[:6])
        )
    typer.echo(f"{structure_path.name} chain A matches P02786 at all {len(residues)} "
               f"modelled residues")

    scan = read_patches(patches)
    measured = read_measured_contacts(contacts)
    thresholds_path = patches.with_suffix(".thresholds.json")
    thresholds = json.loads(thresholds_path.read_text()) if thresholds_path.exists() else {}
    typer.echo(f"{len(scan)} usable patches and {len(measured)} measured contacts read\n")

    centres = {r.seqid.num: c for r in residues if (c := centroid(r)) is not None}
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = {
        name: build_entry(
            name, epitope, check_epitope_against_scan(name, epitope, scan, measured),
            residues, centres, sequence, out_dir, measured,
        )
        for name, epitope in EPITOPES.items()
    }

    manifest = {
        "source": structure_path.name,
        "cut_radius": CUT_RADIUS,
        "min_run": MIN_RUN,
        "acidic_residues": sorted(ACIDIC),
        "scan_thresholds": thresholds,
        "epitopes": entries,
    }
    (out_dir / "hotspots.json").write_text(json.dumps(manifest, indent=2) + "\n")
    typer.echo(f"\nwrote {out_dir / 'hotspots.json'}")


if __name__ == "__main__":
    app()
