"""Verify a pool of designed metal binders and log the filtering to Weights & Biases.

Demonstrates the verifier's value prop: generate many candidates (RFdiffusionAA →
LigandMPNN packs), let the geometry verifier filter to the trustworthy few. The W&B
wiring lives in `touchstone.tracking`; ingestion in `touchstone.load_designs`.

    uv run --extra viz python scripts/log_wandb.py <design_pdb_dir>

Reads WANDB_API_KEY / WANDB_PROJECT from touchstone/.env (gitignored).
"""

from __future__ import annotations

import typer
from rich.console import Console
from rich.table import Table

from touchstone import GeometryVerifier, PDBReference, load_designs, rank, tracking

console = Console()


def main(design_dir: str = typer.Argument("ligmpnn_out")) -> None:
    ref = PDBReference()
    with console.status("verifying designs…"):
        designs = load_designs(f"{design_dir}/**/*.pdb", generator="rfaa+ligmpnn")
        ranked = rank(designs, GeometryVerifier(ref))

    run = tracking.init("ligmpnn-nickel-filter", config={"n_designs": len(ranked), "metal": "Ni2+"})
    counts = tracking.log_ranked(run, ranked, ref)

    # 3D structures: best overall + the best of each verdict class (d.source = its PDB path)
    mols, seen = {"best_design_3d": ranked[0][0].source}, set()
    for d, v in ranked:
        if v.label not in seen:
            seen.add(v.label)
            mols[f"3d_best_{v.label}"] = d.source
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
