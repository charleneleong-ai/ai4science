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
def test_eval(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    test: Path = typer.Argument(..., help="path to test.parquet"),
    ranker: str = typer.Option("fragments", help=" | ".join(RANKERS)),
    ppm: float = typer.Option(5.0, help="neutral-mass retrieval window"),
) -> None:
    """Score a ranker on the REAL test molecules, without using the duplication.

    98% of test spectra are bit-identical to a training row, which hands over the labels.
    That makes an honest evaluation possible rather than impossible: use the duplication only
    to recover ground truth, then score a ranker that never reads a training spectrum. The
    fragment ranker qualifies — it sees a candidate's structure and the query spectrum, and
    nothing else — so this is the best available estimate of performance on molecules the
    leak did not solve for us.
    """
    from .submit import build_lookup, read_test
    from .harness import reciprocal_rank, retrieve

    console = Console()
    chosen = load_ranker(ranker)
    with console.status("recovering test labels via the duplicate join..."):
        labels = build_lookup(train, want_key=True)
    molecules = read_test(test)
    with console.status("building the structure pool..."):
        structures = build_pool(train)
    masses = [s.mass for s in structures]

    scored = hits = unresolved = 0
    rr_total = 0.0
    with console.status(f"scoring {len(molecules)} test molecules with {ranker!r}..."):
        for mol in molecules:
            truth = next((labels[k] for k in mol.keys if k in labels), None)
            if truth is None:
                unresolved += 1
                continue
            query = mol.as_query()
            if query.mass is None:
                continue
            scored += 1
            candidates = retrieve(structures, masses, query.mass, ppm)
            if any(s.inchikey14 == truth for s in candidates):
                hits += 1
                rr_total += reciprocal_rank(chosen(query, candidates), truth)

    console.print()
    console.print(f"ranker               {ranker}")
    console.print(f"test molecules       {len(molecules)}  (labels unrecovered: {unresolved})")
    console.print(f"scored               {scored}")
    console.print(f"candidate recall     {hits / max(scored, 1):.4f}   <- ceiling")
    console.print(f"MRR@25               {rr_total / max(scored, 1):.4f}   <- leak-free estimate")


@app.command()
def submit(
    train: Path = typer.Argument(..., help="path to train.parquet"),
    test: Path = typer.Argument(..., help="path to test.parquet"),
    out: Path = typer.Option(Path("submission.csv"), help="where to write the submission"),
    ranker: str = typer.Option("fragments", help="ranker for molecules the lookup misses"),
    ppm: float = typer.Option(5.0, help="neutral-mass retrieval window for the fallback"),
    use_lookup: bool = typer.Option(
        True, "--lookup/--no-lookup",
        help="--no-lookup ignores the duplication and submits only what the ranker predicts",
    ),
) -> None:
    """Write a submission, reporting how much came from the duplication versus from ranking.

    The split matters more than the score. With `--lookup` (default) the join answers every
    molecule, so the score measures a join rather than a model. With `--no-lookup` the
    submission is what the ranker actually predicts — measured at MRR@25 0.66 on these same
    molecules, so honesty is not expensive here.
    """
    from .submit import build_lookup, pad, rank_fallback, read_test, write_submission

    console = Console()
    lookup: dict = {}
    if use_lookup:
        with console.status("indexing train for the duplicate join..."):
            lookup = build_lookup(train)
    molecules = read_test(test)
    console.print(f"test molecules: [bold]{len(molecules)}[/]   train join keys: {len(lookup):,}")

    rows: dict[str, list[str]] = {}
    from_lookup = 0
    needs_ranking: list = []
    for mol in molecules:
        hits = [lookup[k] for k in mol.keys if k in lookup]
        if hits:
            from_lookup += 1
            # dedupe while keeping order; a molecule's spectra may resolve to one structure
            rows[mol.molecule_id] = pad(list(dict.fromkeys(hits)))
        else:
            needs_ranking.append(mol)

    if needs_ranking:
        chosen = load_ranker(ranker)
        with console.status(f"ranking {len(needs_ranking)} unresolved molecules..."):
            structures = build_pool(train)
            masses = [s.mass for s in structures]
            for mol in needs_ranking:
                rows[mol.molecule_id] = pad(rank_fallback(mol, structures, masses, chosen, ppm))

    write_submission(rows, out)
    console.print()
    console.print(f"from the duplicate join  [bold]{from_lookup}[/] / {len(molecules)}")
    console.print(f"from the {ranker} ranker  {len(needs_ranking)} / {len(molecules)}")
    console.print(f"wrote {out}  ({len(rows)} rows x {len(next(iter(rows.values())))} candidates)")
    if from_lookup:
        console.print("[yellow]note:[/] the join reflects test spectra duplicated in train, "
                      "so any score it earns measures a lookup rather than a model.")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
