"""Molecule-disjoint evaluation, from the command line."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from .harness import (
    Ranker,
    Structure,
    build_coconut_pool,
    build_pool,
    collect_queries,
    evaluate,
    merge_pools,
    rank_by_mass_error,
)

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


def assemble_pool(kind: str, train: Path, coconut: Path | None) -> list[Structure]:
    """`train` holds the answers because they share a file; `coconut` is the realism check."""
    if kind == "train":
        return build_pool(train)
    if coconut is None:
        raise typer.BadParameter(f"pool={kind!r} needs --coconut pointing at the COCONUT csv")
    external = build_coconut_pool(coconut)
    if kind == "coconut":
        return external
    if kind == "both":
        return merge_pools(build_pool(train), external)
    raise typer.BadParameter(f"unknown pool {kind!r}; choose train, coconut or both")


@app.command()
def split_eval(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    ranker: str = typer.Option("mass", help=" | ".join(f"{k}: {v}" for k, v in RANKERS.items())),
    pool: str = typer.Option("train", help="train (holds the answers) | coconut (realism) | both"),
    coconut: Path = typer.Option(None, help="COCONUT csv, required for pool=coconut|both"),
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
    with console.status(f"building the {pool!r} structure pool..."):
        structures = assemble_pool(pool, train, coconut)
    console.print(f"pool ({pool}): [bold]{len(structures):,}[/] unique structures")

    with console.status("collecting held-out queries..."):
        queries = collect_queries(train, frac, limit)
    console.print(f"held out: [bold]{len(queries):,}[/] molecules "
                  f"({sum(len(q.spectra) for q in queries):,} spectra)")

    with console.status(f"scoring with the {ranker!r} ranker..."):
        report = evaluate(queries, structures, ppm=ppm, ranker=chosen)
    console.print()
    console.print(f"pool                 {pool} ({len(structures):,})")
    console.print(f"ranker               {ranker}")
    console.print(report.render())


@app.command()
def submit(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    test: Path = typer.Argument(..., help="path to test.parquet"),
    out: Path = typer.Option(Path("submission.csv"), help="where to write the submission"),
    ranker: str = typer.Option("fragments", help=" | ".join(RANKERS)),
    ppm: float = typer.Option(5.0, help="neutral-mass retrieval window"),
) -> None:
    """Write a local submission CSV from the ranker over a train-structure pool.

    Local only: CASMI is a code competition, so this file cannot be uploaded — the Kaggle
    notebook in notebooks/ is the submission path. And the pool is train's own structures,
    which do not contain most of the competition's answers: this configuration scored 0.115
    on the leaderboard. A realistic submission needs a PubChem-scale pool.
    """
    from .submit import pad, rank_fallback, read_test, write_submission

    console = Console()
    chosen = load_ranker(ranker)
    molecules = read_test(test)
    with console.status("building the structure pool..."):
        structures = build_pool(train)
        masses = [s.mass for s in structures]
    with console.status(f"ranking {len(molecules)} molecules with {ranker!r}..."):
        rows = {
            mol.molecule_id: pad(rank_fallback(mol, structures, masses, chosen, ppm))
            for mol in molecules
        }
    write_submission(rows, out)
    console.print(f"wrote {out}  ({len(rows)} rows x {len(next(iter(rows.values())))} candidates)")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
