"""Build the single ranked submission CSV in Proteinbase's actual round-2 format.

Takes the per-track shortlists `assemble_submission.py` writes and merges them into the one
file that is uploaded. The format rules -- accepted molecule classes, the 20-row cap, the
10-250 aa bounds -- are `tnf_binder.submission`'s, which is also what `tnf-binder validate`
checks the result against, so a file cannot pass the writer and fail the checker.

Ranking is a single combined order, top row first, since there is one file rather than one per
track: pH switch first (the top-ranked objective, and the FAQ says a weak-but-selective binder
beats a strong non-selective one), then interface confidence, then cross-reactivity.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.submission import (MAX_ROWS, MOLECULE_CLASS, SUBMISSION_COLUMNS, row_problem,
                                   write_submission)
from tnf_binder.tables import number, read_rows

app = typer.Typer(add_completion=False)
console = Console()


@app.command()
def main(
    tracks: list[Path] = typer.Argument(..., help="per-track CSVs to merge"),
    out: Path = typer.Option(..., "--out", help="the file that gets uploaded"),
    cap: int = typer.Option(MAX_ROWS, "--cap", help="rows per submission"),
) -> None:
    """Merge, rank and cap the shortlists into one upload-ready CSV."""
    seen, rows = set(), []
    for track in tracks:
        for row in read_rows(track):
            sequence = row["sequence"].strip().upper()
            if sequence in seen:
                continue
            seen.add(sequence)
            rows.append(row)

    rows.sort(key=lambda r: (-number(r, "ph_switch_score"), -number(r, "boltz_iptm_mean"),
                             -number(r, "cross_reactivity")))
    rows = rows[:cap]

    named = [(f"tnf_{i:02d}", r["sequence"].strip().upper()) for i, r in enumerate(rows, 1)]
    for name, sequence in named:
        if problem := row_problem({"name": name, "sequence": sequence,
                                   "molecule_class": MOLECULE_CLASS}):
            typer.echo(f"refusing to write: {problem}", err=True)
            raise typer.Exit(code=1)

    out.parent.mkdir(parents=True, exist_ok=True)
    write_submission(out, SUBMISSION_COLUMNS,
                     [[name, sequence, MOLECULE_CLASS] for name, sequence in named])
    console.print(f"wrote [bold]{len(named)}[/] designs -> [cyan]{out}[/] (cap {cap})")
    for (name, sequence), row in zip(named, rows):
        console.print(f"  {name} len={len(sequence):3d} ph={row['ph_switch_score']:>4} "
                      f"iptm={row['boltz_iptm_mean']:>6} xr={row['cross_reactivity']:>7} "
                      f"his={row['pka_detail'] or '-'}")


if __name__ == "__main__":
    app()
