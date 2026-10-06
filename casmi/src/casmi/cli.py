"""Molecule-disjoint evaluation, from the command line."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from .harness import Ranker, build_pool, collect_queries, evaluate, rank_by_mass_error

app = typer.Typer(add_completion=False, help="CASMI 2026 molecule-disjoint harness.")

RANKERS: dict[str, str] = {
    "mass": "|mass error| only — the floor a real ranker must beat",
    "fragments": "combinatorial fragmentation: share of peak intensity the structure explains",
}


def load_ranker(name: str) -> Ranker:
    """rdkit is only imported for the ranker that needs it, so `mass` runs without it."""
    if name == "mass":
        return rank_by_mass_error
    if name == "fragments":
        from .fragments import rank_by_fragments

        return rank_by_fragments
    raise typer.BadParameter(f"unknown ranker {name!r}; choose from {sorted(RANKERS)}")


@app.command()
def split_eval(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    ranker: str = typer.Option("mass", help=" | ".join(f"{k}: {v}" for k, v in RANKERS.items())),
    frac: float = typer.Option(0.02, help="fraction of molecules held out (hashed on inchikey14, so stable)"),
    limit: int = typer.Option(500, help="cap on held-out molecules scored (0 = all)"),
    ppm: float = typer.Option(5.0, help="neutral-mass retrieval window"),
) -> None:
    """Measure candidate recall and the MRR@25 floor on held-out molecules.

    Recall and ranking are reported separately: if recall is the binding constraint, no
    amount of ranking work can help, and that is worth knowing before building a ranker.
    """
    console = Console()
    chosen = load_ranker(ranker)
    with console.status("building the structure pool..."):
        pool = build_pool(train)
    console.print(f"pool: [bold]{len(pool):,}[/] unique structures")

    with console.status("collecting held-out queries..."):
        queries = collect_queries(train, frac, limit)
    console.print(f"held out: [bold]{len(queries):,}[/] molecules "
                  f"({sum(len(q.spectra) for q in queries):,} spectra)")

    with console.status(f"scoring with the {ranker!r} ranker..."):
        report = evaluate(queries, pool, ppm=ppm, ranker=chosen)
    console.print()
    console.print(f"ranker               {ranker}")
    console.print(report.render())


def main() -> None:
    app()


if __name__ == "__main__":
    main()
