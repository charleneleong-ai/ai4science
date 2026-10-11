"""Collecting the generator's own confidence, and joining it to what we measured ourselves.

    metrics.py collect <fold_out/> --out folding.csv
    metrics.py join folding.csv geometry.csv --out merged.csv

Kept as two stages because the generator's score and ours are different evidence and must stay
separable. BoltzGen grading its own fold is not independent of the design it produced, so it is
collected, reported, and then joined to a geometric score computed from the structure -- never
blended into one number. Where they disagree is the interesting part: AlphaFold2 rated its own
optimised designs at i_pTM 0.92 where an independent predictor rated the same structures 0.21.

The join also gates, and gating on both sources is the point: a design can fold confidently and
still bury almost nothing on the epitope.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import typer
from rich.console import Console

from tnf_binder.tables import number, read_rows, write_rows

app = typer.Typer(add_completion=False)
console = Console()


# stdout is the CSV itself; notes belong on stderr so the table can be redirected to a file.

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
def collect(
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


COLUMNS = ["design", "binder_len", "iptm", "design_to_target_iptm", "interaction_pae",
           "design_ptm", "target_ptm", "interface_area", "epitope_focus", "cross_reactivity",
           "divergent_frac", "clashes", "flag", "sequence"]
NO_PAE = 99.0


def merged_row(scored: dict[str, str], folded: dict[str, str], min_d2t: float,
               max_pae: float) -> dict[str, object]:
    d2t = number(folded, "design_to_target_iptm")
    pae = number(folded, "interaction_pae", NO_PAE)
    flags = [f for f in (scored.get("flag") or "").split(",") if f]
    if not folded:
        flags.append("no_fold_metrics")
    else:
        if d2t < min_d2t:
            flags.append(f"low_d2t_iptm({d2t:.2f})")
        if pae > max_pae:
            flags.append(f"high_pae({pae:.1f})")
    return scored | {
        "iptm": number(folded, "iptm"), "design_to_target_iptm": d2t, "interaction_pae": pae,
        "target_ptm": number(folded, "target_ptm"), "design_ptm": number(folded, "design_ptm"),
        "flag": ",".join(flags),
    }


@app.command()
def join(
    scored_csv: Path = typer.Argument(..., help="geometry table from score_designs.py"),
    fold_csv: Path = typer.Argument(..., help="confidence table from read_fold_metrics.py"),
    out: Path = typer.Option(..., "--out", help="the joined, gated table"),
    min_d2t_iptm: float = typer.Option(0.5, "--min-d2t-iptm",
                                       help="threshold on design_to_target_iptm, the "
                                            "binder-specific interface score"),
    max_pae: float = typer.Option(10.0, "--max-pae"),
    relax: bool = typer.Option(False, "--relax",
                               help="keep designs failing the confidence gate, flagged, "
                                    "for inspection"),
) -> None:
    """Join the two tables and keep only what clears both gates."""
    folded = {r["design"]: r for r in read_rows(fold_csv)}
    merged = [merged_row(r, folded.get(r["design"], {}), min_d2t_iptm, max_pae)
              for r in read_rows(scored_csv)]

    passing = [r for r in merged if not r["flag"]]
    kept = merged if relax else passing
    kept.sort(key=lambda r: (-number(r, "design_to_target_iptm"),
                             number(r, "interaction_pae", NO_PAE)))
    write_rows(out, kept, COLUMNS)

    console.print(f"{len(merged)} scored | [bold]{len(passing)}[/] pass both gates "
                  f"(design_to_target_iptm>={min_d2t_iptm}, pae<={max_pae}) | "
                  f"wrote {len(kept)} -> [cyan]{out}[/]")
    if not passing and merged:
        best = max(merged, key=lambda r: number(r, "design_to_target_iptm"))
        console.print(f"  best by design_to_target_iptm: {best['design']} "
                      f"d2t={best['design_to_target_iptm']:.3f} "
                      f"(global iptm={best['iptm']:.3f}) "
                      f"pae={best['interaction_pae']:.1f} -- nothing clears the bar")


if __name__ == "__main__":
    app()
