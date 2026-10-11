"""Collect BoltzGen folding-step confidence metrics.

Each metric is stored per DIFFUSION SAMPLE, not as a scalar: coords has shape (5, n_atoms, 3)
and every confidence key is shape (5,). The samples are not interchangeable -- they trend
across the schedule, e.g. design_to_target_iptm [0.237, 0.287, 0.325, 0.352, 0.309] for one
design -- so collapsing them to a single number silently picks a statistic.

Both are reported:
  *_mean  describes the config, for comparing settings
  *_best  the sample you would actually keep, so it is what ranking uses

An earlier version returned the mean for size>1 and the first element elsewhere, which made two
readers disagree on the same file. Reporting both removes the ambiguity.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import typer
from rich.console import Console

from tnf_binder.tables import write_rows

app = typer.Typer(add_completion=False)
# stdout is the CSV itself; notes belong on stderr so the table can be redirected to a file.
console = Console(stderr=True)

# interface confidence: higher is better / error: lower is better
HIGHER = ("iptm", "design_iptm", "design_to_target_iptm", "design_iiptm",
          "ptm", "design_ptm", "target_ptm", "protein_iptm")
LOWER = ("interaction_pae", "min_interaction_pae", "min_design_to_target_pae")


def metrics_of(design: str, archive: np.lib.npyio.NpzFile) -> dict[str, object]:
    """One row per npz, with both statistics for every confidence key it carries."""
    row: dict[str, object] = {"design": design}
    for key in HIGHER + LOWER:
        if key not in archive.files:
            continue
        values = np.asarray(archive[key]).reshape(-1).astype(float)
        if not values.size:
            continue
        row[f"{key}_mean"] = round(float(values.mean()), 4)
        row[f"{key}_best"] = round(
            float(values.max() if key in HIGHER else values.min()), 4)
    return row


def sample_count(archive: np.lib.npyio.NpzFile) -> int | None:
    if "coords" not in archive.files:
        return None
    return int(np.asarray(archive["coords"]).shape[0])


@app.command()
def main(
    npz_dir: Path = typer.Argument(..., help="a BoltzGen fold_out_npz directory"),
    out: Path | None = typer.Option(None, "--out", help="also write the table here"),
    sort: str = typer.Option("design_to_target_iptm_best", "--sort",
                             help="column to rank on, highest first"),
) -> None:
    """Read every npz under a directory and report mean and best per confidence key."""
    paths = sorted(npz_dir.rglob("*.npz"))
    if not paths:
        typer.echo(f"no npz under {npz_dir}", err=True)
        raise typer.Exit(code=2)

    rows, samples = [], set()
    for path in paths:
        try:
            archive = np.load(path, allow_pickle=True)
            rows.append(metrics_of(path.stem, archive))
            samples.add(sample_count(archive))
        except Exception as exc:  # one unreadable archive must not sink the batch
            console.print(f"[dim]skip {path.name}: {type(exc).__name__}[/]")
    if counts := sorted(n for n in samples if n is not None):
        console.print(f"[dim]diffusion samples per design: {counts}[/]")

    columns = ["design"] + [c for c in
                            [f"{k}_{s}" for k in HIGHER + LOWER for s in ("best", "mean")]
                            if any(c in r for r in rows)]
    rows.sort(key=lambda r: -float(r.get(sort) or -1e9))
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    if out:
        write_rows(out, rows, columns)
        console.print(f"\nwrote [cyan]{out}[/]")


if __name__ == "__main__":
    app()
