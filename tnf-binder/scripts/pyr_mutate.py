"""Mutate one residue to histidine with a real, repacked side chain, via PyRosetta.

A single-site front end onto `tnf_binder.rosetta.HisMutator`, kept because the His scans are
long-running and a one-off mutation is how a failing site gets inspected by hand.

Runs in the BindCraft conda env, which is where pyrosetta lives. PYTHONNOUSERSITE=1 is required:
~/.local numpy shadows that env otherwise.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.rosetta import HisMutator

app = typer.Typer(add_completion=False)
console = Console()


@app.command()
def main(
    pdb: Path = typer.Argument(..., help="the complex to mutate"),
    chain: str = typer.Option(..., "--chain"),
    resnum: int = typer.Option(..., "--resnum"),
    out: Path = typer.Option(..., "--out"),
    repack_radius: float = typer.Option(6.0, "--repack-radius"),
) -> None:
    """Substitute one position with histidine and write the repacked structure."""
    mutator = HisMutator(repack_radius=repack_radius)
    pose = mutator.load(pdb)
    if not mutator.mutate(pose, [(chain, resnum)]):
        typer.echo(f"{chain}{resnum} is not in the pose, or did not become HIS", err=True)
        raise typer.Exit(code=2)
    mutator.write(pose, out)
    console.print(f"[green]OK[/] {chain}{resnum} -> HIS; wrote [cyan]{out}[/]")


if __name__ == "__main__":
    app()
