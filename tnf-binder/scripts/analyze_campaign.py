"""Does histidine engineering help binding? Paired parent-vs-variant test across designs.

On one design all four variants beat the parent, but with high per-sample variance
(0.586 / 0.655 / 0.192 for one variant). This pairs every variant against its OWN parent,
folded under identical conditions, so design-to-design differences cancel and the question
becomes a within-family comparison.

Reports the binder-to-target ipTM only -- global ipTM is inflated by the TNF trimer's own
internal interfaces (~0.93) and is useless for ranking binders. The mean across diffusion
samples is the statistic; the best sample is shown beside it for context only.
"""

from __future__ import annotations

import statistics as stats
from collections import defaultdict
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from tnf_binder.confidence import DesignConfidence
from tnf_binder.tables import read_rows, write_rows

app = typer.Typer(add_completion=False)
console = Console()

COLUMNS = ["spec", "design", "kind", "positions", "pkas", "best", "mean", "ipde", "sequence"]


def measured_rows(manifest: list[dict[str, str]], pred_root: Path) -> list[dict[str, object]]:
    """Manifest rows that have a folded prediction, with their confidence attached."""
    rows: list[dict[str, object]] = []
    for row in manifest:
        folded = DesignConfidence.from_dir(pred_root / row["spec"])
        if not folded:
            continue
        rows.append(row | {"best": folded.binder_to_target_best,
                           "mean": folded.binder_to_target,
                           "ipde": folded.complex_ipde})
    return rows


def families(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    """design -> {parent, variants}, so every variant is compared against its own parent."""
    found: dict[str, dict[str, object]] = defaultdict(dict)
    for row in rows:
        if row["kind"] == "parent":
            found[row["design"]]["parent"] = row
        else:
            found[row["design"]].setdefault("variants", []).append(row)
    return found


def paired_table(grouped: dict[str, dict[str, object]]) -> tuple[Table, list[float]]:
    table = Table(title="parent vs its best His variant · mean binder-to-target ipTM")
    for column in ("design", "parent best", "parent mean", "variant best", "variant mean",
                   "delta mean", "best variant"):
        table.add_column(column, justify="right" if "design" not in column else "left")

    deltas = []
    for design, family in sorted(grouped.items()):
        parent, variants = family.get("parent"), family.get("variants", [])
        if not parent or not variants:
            continue
        best = max(variants, key=lambda r: r["mean"])
        delta = best["mean"] - parent["mean"]
        deltas.append(delta)
        colour = "green" if delta > 0 else "red"
        table.add_row(design, f"{parent['best']:.3f}", f"{parent['mean']:.3f}",
                      f"{best['best']:.3f}", f"{best['mean']:.3f}",
                      f"[{colour}]{delta:+.3f}[/]", str(best["positions"]))
    return table, deltas


@app.command()
def main(
    manifest: Path = typer.Argument(..., help="the campaign manifest from build_his_campaign"),
    pred_root: Path = typer.Argument(..., help="Boltz predictions/ root, one dir per spec"),
    out: Path | None = typer.Option(None, "--out", help="write the measured rows here"),
) -> None:
    """Pair every His variant against its own parent and report the sign test."""
    rows = measured_rows(read_rows(manifest), pred_root)
    if not rows:
        typer.echo(f"nothing folded under {pred_root}", err=True)
        raise typer.Exit(code=2)

    table, deltas = paired_table(families(rows))
    console.print(table)
    if deltas:
        improved = sum(1 for d in deltas if d > 0)
        console.print(f"\n[bold]PAIRED RESULT across {len(deltas)} families[/]: "
                      f"variants improved the parent mean in {improved}/{len(deltas)}")
        spread = f", stdev {stats.stdev(deltas):.3f}" if len(deltas) > 1 else ""
        console.print(f"  mean delta {stats.mean(deltas):+.3f}, "
                      f"median {stats.median(deltas):+.3f}{spread}")
        console.print(f"  sign test: {improved} up / {len(deltas) - improved} down")

    console.print("\n[bold]Top 12 constructs by mean binder-to-target ipTM[/]")
    for row in sorted(rows, key=lambda r: -r["mean"])[:12]:
        console.print(f"  {row['spec']:46s} {row['kind']:7s} best={row['best']:.3f} "
                      f"mean={row['mean']:.3f} ipde={row['ipde']:.2f}  {row['pkas']}")

    if out:
        write_rows(out, sorted(rows, key=lambda r: -r["mean"]), COLUMNS)
        console.print(f"\nwrote [cyan]{out}[/]")


if __name__ == "__main__":
    app()
