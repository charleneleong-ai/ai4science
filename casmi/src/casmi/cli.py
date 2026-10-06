"""Molecule-disjoint evaluation, from the command line."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from .harness import build_pool, collect_queries, evaluate

app = typer.Typer(add_completion=False, help="CASMI 2026 molecule-disjoint harness.")


@app.command()
def split_eval(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    frac: float = typer.Option(0.02, help="fraction of molecules held out (hashed on inchikey14, so stable)"),
    limit: int = typer.Option(500, help="cap on held-out molecules scored (0 = all)"),
    ppm: float = typer.Option(5.0, help="neutral-mass retrieval window"),
) -> None:
    """Measure candidate recall and the MRR@25 floor on held-out molecules.

    Recall and ranking are reported separately: if recall is the binding constraint, no
    amount of ranking work can help, and that is worth knowing before building a ranker.
    """
    console = Console()
    with console.status("building the structure pool..."):
        pool = build_pool(train)
    console.print(f"pool: [bold]{len(pool):,}[/] unique structures")

    with console.status("collecting held-out queries..."):
        queries = collect_queries(train, frac, limit)
    console.print(f"held out: [bold]{len(queries):,}[/] molecules "
                  f"({sum(len(q.spectra) for q in queries):,} spectra)")

    report = evaluate(queries, pool, ppm=ppm)
    console.print()
    console.print(report.render())


def main() -> None:
    app()


if __name__ == "__main__":
    main()
