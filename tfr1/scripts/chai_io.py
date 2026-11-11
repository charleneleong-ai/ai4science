"""Both directions of the Chai-1 boundary: write its inputs, then collect its poses.

Why a second model at all is argued in `glp1r/scripts/chai_independent.py`, the driver this
feeds; in short, BoltzGen generates *and* ranks with Boltz-2 weights, so its own interface score
is one model family grading itself.

Collecting the poses back is a conversion rather than a copy because of a layout mismatch, and
the conversion is the shallower of the two available fixes. Chai writes
`pred.model_trunk_0_idx_<k>.pdb` inside a directory named for the job, while
`ph_switch.index_structures` globs `*.cif*` and looks for the design name in the *file stem*. A
flat `<design>_idx<k>.cif` satisfies it. The root fix is in `index_structures`: accept `.pdb`
(`gemmi.read_structure` reads it natively, so only discovery objects) and match the design name
against the path relative to the root, after which a directory-per-design layout indexes in
place and this command deletes. That is deferred because the pose measurement it would serve
turned out to be invalid -- Chai does not reproduce the designed interface on these targets --
so it would be infrastructure for a measurement nobody can use.

Chai ignores the `name=` fields in the FASTA headers and labels the chains in input order, so
the binder is written first.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Annotated

import gemmi
import typer
from csvio import pool_rows, write_rows
from target_sites import chain_sequence

BINDER_CHAIN, TARGET_CHAIN = "X", "R"

app = typer.Typer(add_completion=False)


def chain_sequences(path: Path, *chain_ids: str) -> dict[str, str]:
    """One parse of a structure file, however many of its chains are wanted."""
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    found = {chain.name: chain_sequence(chain) for chain in structure[0]}
    if missing := [name for name in chain_ids if name not in found]:
        raise typer.BadParameter(f"{path} has no chain(s) {missing}: {sorted(found)}")
    return {name: found[name] for name in chain_ids}


def read_chain(path: Path, chain_id: str) -> str:
    return chain_sequences(path, chain_id)[chain_id]


def write_fasta(path: Path, binder: str, target: str) -> None:
    path.write_text(f">protein|name={BINDER_CHAIN}\n{binder}\n"
                    f">protein|name={TARGET_CHAIN}\n{target}\n")


@app.command()
def inputs(
    candidates: Annotated[list[Path], typer.Option(
        help="selection CSV; repeat the option to fold the union of several")],
    targets: Annotated[Path, typer.Option(help="directory holding <epitope>.pdb")],
    out_dir: Annotated[Path, typer.Option(help="where the per-design FASTAs are written")],
    manifest: Annotated[Path, typer.Option(help="what was written, for the measuring step")],
) -> None:
    chosen = pool_rows(candidates)
    sequences = {epitope: read_chain(targets / f"{epitope}.pdb", TARGET_CHAIN)
                 for epitope in sorted({row["epitope"] for row in chosen.values()})}
    for name, sequence in sequences.items():
        typer.echo(f"target {name:14s} chain {TARGET_CHAIN}  {len(sequence)} residues")

    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, row in sorted(chosen.items()):
        binder = row["sequence"].strip().upper()
        target = sequences[row["epitope"]]
        write_fasta(out_dir / f"{name}.fasta", binder, target)
        rows.append({"name": name, "sequence": binder,
                     "parent": row["parent"], "epitope": row["epitope"],
                     "binder_length": len(binder), "target_length": len(target),
                     "substitutions": row.get("substitutions", ""),
                     "switch_log_units_single_pose": row.get("switch_log_units", "")})

    write_rows(manifest, rows)
    typer.echo(f"\nwrote {len(rows)} FASTAs to {out_dir} and the manifest to {manifest}")


def pose_files(job_dir: Path) -> Iterator[tuple[int, Path]]:
    """Chai's per-sample predictions, each with the sample index it was written under."""
    indexed = ((int(path.stem.rsplit("_", 1)[-1]), path)
               for path in job_dir.glob("pred.model_*.pdb"))
    yield from sorted(indexed)


@app.command()
def poses(
    chai_out: Annotated[Path, typer.Option(help="Chai output root, one directory per design")],
    out_dir: Annotated[Path, typer.Option(help="flat directory of <design>_idx<k>.cif")],
) -> None:
    """Flatten Chai's per-job directories into the layout the switch verifier indexes."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    for job_dir in sorted(p for p in chai_out.iterdir() if p.is_dir()):
        samples = list(pose_files(job_dir))
        if not samples:
            typer.echo(f"{job_dir.name:56s} NO SAMPLES - nothing to measure")
            continue
        for index, source in samples:
            structure = gemmi.read_structure(str(source))
            structure.setup_entities()
            structure.make_mmcif_document().write_file(
                str(out_dir / f"{job_dir.name}_idx{index}.cif"))
        written += len(samples)
        typer.echo(f"{job_dir.name:56s} {len(samples)} samples")

    typer.echo(f"\nwrote {written} poses to {out_dir}")


@app.command()
def control(
    crystal: Annotated[Path, typer.Option(help="6WRV, which holds three copies of the complex")],
    cut_target: Annotated[Path, typer.Option(help="the carved target the designs were made for")],
    out_dir: Annotated[Path, typer.Option(help="where the control FASTAs are written")],
    binder_chain: Annotated[str, typer.Option(help="3DS18's chain in the crystal")] = "C",
    native_chain: Annotated[str, typer.Option(help="its TfR1 partner in the crystal")] = "B",
) -> None:
    """A known 6WRV binder against both the native and the carved target.

    Chai scoring the designs near 0.1 admits two readings -- the designs do not bind, or a target
    carved to a surface patch is not something Chai can score. A crystallographic binder for this
    receptor separates them, because its answer is already known. Neither arm is evidence about a
    design; they calibrate the instrument that judges them.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    chains = chain_sequences(crystal, binder_chain, native_chain)
    binder, native = chains[binder_chain], chains[native_chain]
    cut = read_chain(cut_target, TARGET_CHAIN)

    for name, target in (("3ds18-native", native), ("3ds18-cut", cut)):
        write_fasta(out_dir / f"{name}.fasta", binder, target)
        typer.echo(f"{name:16s} binder {len(binder)} res + target {len(target)} res")
    typer.echo(f"\nwrote 2 control FASTAs to {out_dir}")


if __name__ == "__main__":
    app()
