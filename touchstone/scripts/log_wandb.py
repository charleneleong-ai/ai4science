"""Verify a pool of designed metal binders and log the filtering to Weights & Biases.

Demonstrates the verifier's value prop: generate many candidates (RFdiffusionAA →
LigandMPNN packs), let the geometry verifier filter to the trustworthy few. The W&B
wiring lives in `touchstone.tracking`; this script just loads the designs, ranks them,
and hands them over.

    uv run --extra viz python scripts/log_wandb.py <design_pdb_dir>

Reads WANDB_API_KEY / WANDB_PROJECT from touchstone/.env (gitignored).
"""

from __future__ import annotations

import glob
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import track as progress
from rich.table import Table

from touchstone import (
    BinderDesign,
    GeometryVerifier,
    PDBReference,
    coordination_site_from_pdb,
    rank,
    tracking,
)

console = Console()


def main(design_dir: str = typer.Argument("ligmpnn_out")) -> None:
    ref = PDBReference()
    designs = []
    for pdb in progress(sorted(glob.glob(f"{design_dir}/**/*.pdb", recursive=True)),
                        description="verifying designs", console=console):
        try:
            site = coordination_site_from_pdb(pdb, "NI", "Ni2+")
        except ValueError:
            continue
        d = BinderDesign(Path(pdb).stem, site, "rfaa+ligmpnn", float("nan"))
        d.path = pdb  # remember source for the 3D view
        designs.append(d)
    ranked = rank(designs, GeometryVerifier(ref))

    run = tracking.init("ligmpnn-nickel-filter", config={"n_designs": len(ranked), "metal": "Ni2+"})
    counts = tracking.log_ranked(run, ranked, ref)

    # 3D structures: best overall + the best of each verdict class
    mols, seen = {"best_design_3d": ranked[0][0].path}, set()
    for d, v in ranked:
        cls = "defer" if v.ood else ("trust" if v.trust else "weak")
        if cls not in seen:
            seen.add(cls)
            mols[f"3d_best_{cls}"] = d.path
    tracking.log_molecules(run, mols)

    table = Table(title=f"Verifier filtering — {len(ranked)} designs")
    table.add_column("verdict")
    table.add_column("count", justify="right")
    for verdict, style in [("trust", "green"), ("weak", "yellow"), ("defer", "red")]:
        table.add_row(verdict, str(counts[verdict]), style=style)
    console.print(table)
    console.print(f"[bold]logged →[/] {run.url}")
    run.finish()


if __name__ == "__main__":
    typer.run(main)
