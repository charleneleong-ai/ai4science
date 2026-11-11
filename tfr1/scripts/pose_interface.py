"""Where a co-folding model actually puts the binder, which its ipTM does not say.

An interface ipTM near 0.1 is read as "this design does not bind", but it cannot distinguish that
from "the model docked the design somewhere else". The two have opposite consequences: the first
is a verdict on the design, the second is a verdict on the prediction, and only the first is
evidence about a submission. The difference is in the contacts -- which target residues the
binder touches, whether they are the hotspots the generator was given, and whether the model's
own samples agree with each other.

Three things are refused rather than guessed, because each would put a plausible wrong number in
a column that travels with a submission. The binder chain is found by matching the recorded
sequence, as `ph_switch.load_complex` does, because Chai relabels chains and ignores the `name=`
fields in its FASTA input, so neither chain order nor chain name survives the round trip. A
complex with more than one target chain is an error, not a silent measurement against the first
of them. And the target chain's sequence is checked against the carved target the hotspot
positions were written for, since `hotspots_positions` is a 1-based index into that chain
(`prepare_target.verify_positions`) -- one dropped residue shifts every position and nothing
downstream would notice.

Contacts are a naive all-pairs scan at `target_sites.CONTACT_CUTOFF`, the repo's one definition
of a heavy-atom contact. `target_sites.contact_numbers` already does this with
`gemmi.NeighborSearch`, carrying a hard-won note about symmetry images on 6WRV (which does not
arise here -- these poses are P1 and the search would run against an empty unit cell either
way), and the two should share one helper.

That sharing is deferred, and the honest cost of deferring it is ~30 s per full run, not the
"seconds once" an earlier version of this note claimed. Benchmarked on these poses: the scan is
~25 ms each against ~2.6 ms for one `NeighborSearch` answering both questions, so ~33 s against
~3.4 s over a 1300-pose batch, and the gap grows with the pool rather than being a one-off.
Index construction is 0.03 ms, so nothing has to be cached across poses. It is deferred anyway
because the measurement is already made and frozen for this submission, and re-cutting its
inner loop the day before the deadline risks the numbers it produced to save half a minute of
compute -- not because the saving is small.
"""

from __future__ import annotations

import json
from collections.abc import Collection, Iterator, Sequence
from pathlib import Path
from statistics import mean
from typing import Annotated

import gemmi
import typer
from chai_io import TARGET_CHAIN, pose_files
from csvio import pool_rows, same_design, write_rows
from target_sites import CONTACT_CUTOFF, HYDROGEN, chain_sequence, load_chain

app = typer.Typer(add_completion=False)


def binder_and_target(path: Path, sequence: str) -> tuple[gemmi.Chain, gemmi.Chain]:
    model = gemmi.read_structure(str(path))[0]
    binders = [chain for chain in model
               if len(chain) == len(sequence) and chain_sequence(chain) == sequence]
    if len(binders) != 1:
        shapes = {chain.name: len(chain) for chain in model}
        raise typer.BadParameter(
            f"{path.name}: {len(binders)} chains carry the recorded {len(sequence)}-residue "
            f"sequence; chains are {shapes}. Refusing to guess which one is the binder.")
    targets = [chain for chain in model if chain.name != binders[0].name]
    if len(targets) != 1:
        raise typer.BadParameter(
            f"{path.name}: {len(targets)} target chains, so the hotspot positions have no one "
            f"chain to index into.")
    return binders[0], targets[0]


def heavy_atoms(chain: gemmi.Chain) -> list[gemmi.Position]:
    """Every non-hydrogen atom's position, flattened once.

    `gemmi` builds a fresh wrapper on each `atom.pos` access, so walking the hierarchy inside a
    distance loop costs more than the distances; both scans below take this list instead. It is
    rebuilt per pose on purpose -- the five Chai samples of one design differ in coordinates.
    """
    return [atom.pos for residue in chain for atom in residue if atom.element != HYDROGEN]


def contacting_positions(binder_atoms: Sequence[gemmi.Position], target: gemmi.Chain,
                         cutoff: float) -> frozenset[int]:
    """The target's 1-based chain positions within `cutoff` of any binder heavy atom."""
    touched = set()
    for position, residue in enumerate(target, start=1):
        for atom in residue:
            if atom.element == HYDROGEN:
                continue
            # Bound once per atom: `atom.pos` builds a fresh wrapper on each access, and this is
            # the innermost loop of the whole script.
            here = atom.pos
            if any(here.dist(other) <= cutoff for other in binder_atoms):
                touched.add(position)
                break
    return frozenset(touched)


def nearest_hotspot_distance(binder_atoms: Sequence[gemmi.Position], target: gemmi.Chain,
                             hotspots: Collection[int]) -> float:
    """Closest approach, in angstroms, between the binder and any nominated hotspot residue.

    The hotspot *count* cannot tell a binder sitting just off the nominated set from one on the
    far side of the carved patch -- both read 0 -- and those two cases mean opposite things
    about whether the nomination or the design is at fault. This distance separates them.
    """
    return min((atom.pos.dist(other)
                for position, residue in enumerate(target, start=1) if position in hotspots
                for atom in residue if atom.element != HYDROGEN
                for other in binder_atoms),
               default=float("nan"))


def find_design(name: str, recorded: dict[str, dict[str, str]]) -> dict[str, str] | None:
    """The table row for a design, tolerating the two spellings of its name.

    A pose written to disk as `conserved_s_45.cif` is recorded in the pooled table as
    `conserved_s_conserved_s_45`; `csvio.same_design` owns that rule. Ambiguity is an error
    rather than a pick, since the wrong row would supply the wrong sequence and the binder
    would then fail to match at all -- or, worse, match and be measured against another
    design's epitope.
    """
    if (exact := recorded.get(name)) is not None:
        return exact
    # Both directions, because `same_design` only tests the file stem inside the recorded name
    # and the poses here are written as `dn_<design>_refold_cifp1.cif`, which is the other way
    # round -- the direction `ph_switch.index_structures` tests. Swapping the arguments reuses
    # the rule rather than restating its underscore padding for a fourth time.
    matches = [row for key, row in recorded.items()
               if same_design(key, name) or same_design(name, key)]
    if len(matches) > 1:
        raise typer.BadParameter(
            f"{name} matches {len(matches)} rows across the design tables, so its sequence "
            f"is ambiguous.")
    return matches[0] if matches else None


def agreement(faces: Sequence[frozenset[int]]) -> float | None:
    """Mean Jaccard overlap of the contacted residues over every pair of poses, or None for one.

    One binding mode sampled repeatedly scores near 1. None is not 0: a single pose has nothing
    to agree with, and reporting 0 there would read as "the samples disagree".
    """
    overlaps = [len(a & b) / len(a | b)
                for index, a in enumerate(faces) for b in faces[index + 1:] if a | b]
    return mean(overlaps) if overlaps else None


def epitope_for(design: dict[str, str], epitopes: Sequence[str]) -> str | None:
    """The epitope the table recorded, or the longest one its name starts with.

    The variant table (`denovo_designs.csv`) carries no epitope column, so for those the name
    prefix is the only record there is. Derivation is counted and reported rather than assumed,
    because a name that does not carry its epitope would otherwise drop out of the funnel.
    """
    if recorded := design.get("epitope"):
        return recorded
    candidates = [name for name in epitopes if design["name"].startswith(name)]
    return max(candidates, key=len) if candidates else None


def jobs(generated: Path | None,
         chai_out: Path | None) -> Iterator[tuple[str, str, list[Path]]]:
    """Each design to measure, as its name, which model produced it, and its poses."""
    for path in sorted(generated.glob("*.cif")) if generated else []:
        yield path.stem, "generator", [path]
    for job_dir in sorted(p for p in chai_out.iterdir() if p.is_dir()) if chai_out else []:
        yield job_dir.name, "chai", [path for _, path in pose_files(job_dir)]


def carved_targets(epitopes: dict[str, dict[str, object]], directory: Path,
                   chain_id: str) -> dict[str, str]:
    """Each epitope's carved target sequence, which its hotspot positions index into."""
    return {name: chain_sequence(load_chain(directory / str(spec["file"]), chain_id))
            for name, spec in epitopes.items()}


@app.command()
def main(
    hotspots: Annotated[Path, typer.Option(help="hotspots.json, for each epitope's binding set")],
    targets: Annotated[Path, typer.Option(help="directory of the carved <epitope>.pdb targets")],
    designs: Annotated[list[Path], typer.Option(
        help="table with name,sequence,epitope; repeat for the chassis and the variants")],
    out: Annotated[Path, typer.Option(help="one row per design per source")],
    chai_out: Annotated[Path | None, typer.Option(
        help="Chai output root, one directory of pred.model_*.pdb per design")] = None,
    generated: Annotated[Path | None, typer.Option(
        help="directory of the generator's own complexes, named <design>.cif")] = None,
    target_chain: Annotated[str, typer.Option(help="the carved target's chain")] = TARGET_CHAIN,
    cutoff: Annotated[float, typer.Option(help="heavy-atom contact distance, angstroms")]
    = CONTACT_CUTOFF,
) -> None:
    epitopes = json.loads(hotspots.read_text())["epitopes"]
    sites = {name: set(spec["hotspots_positions"]) for name, spec in epitopes.items()}
    written = carved_targets(epitopes, targets, target_chain)
    # Epitope from the table that recorded it, never re-derived from the name: a name that does
    # not carry its epitope as a prefix would otherwise drop out of the funnel silently.
    recorded = pool_rows(designs)

    rows: list[dict[str, object]] = []
    derived = 0
    known = list(sites)
    for name, source, poses in jobs(generated, chai_out):
        design = find_design(name, recorded)
        epitope = epitope_for(design, known) if design else None
        if design is None or epitope is None or not poses:
            typer.echo(f"{name:56.56s} {source:10s} SKIPPED ("
                       f"{'no poses' if epitope else 'no recorded sequence or epitope'})")
            continue
        derived += not design.get("epitope")
        sequence = design["sequence"].strip().upper()
        hotspot_set = sites[epitope]
        faces, approaches = [], []
        for path in poses:
            binder, target = binder_and_target(path, sequence)
            if chain_sequence(target) != written[epitope]:
                raise typer.BadParameter(
                    f"{path.name}: its target chain is not the carved {epitope} target the "
                    f"hotspot positions index into, so the join would be off by an unknown "
                    f"shift.")
            atoms = heavy_atoms(binder)
            faces.append(contacting_positions(atoms, target, cutoff))
            approaches.append(nearest_hotspot_distance(atoms, target, hotspot_set))
        overlap = agreement(faces)
        row = {"design": design["name"], "name": name,
               "epitope": epitope, "source": source,
               "poses": len(faces),
               "interface_residues": round(mean(len(f) for f in faces), 1),
               "hotspots_touched": round(mean(len(f & hotspot_set) for f in faces), 2),
               "hotspots_total": len(hotspot_set),
               "nearest_hotspot_angstroms": round(mean(approaches), 1),
               "pose_agreement": "" if overlap is None else round(overlap, 3)}
        rows.append(row)
        typer.echo(f"{row['design']:56.56s} {source:10s} "
                   f"{row['interface_residues']:6} residues  "
                   f"hotspots {row['hotspots_touched']:5}/{row['hotspots_total']} at "
                   f"{row['nearest_hotspot_angstroms']:5}A  "
                   f"agreement {row['pose_agreement']}")

    write_rows(out, rows)
    typer.echo(f"\nwrote {len(rows)} rows to {out}"
               + (f"; {derived} epitopes derived from the design name" if derived else ""))


if __name__ == "__main__":
    app()
