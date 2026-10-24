"""Does a second model agree with the one that made the designs?

BoltzGen generates and ranks with Boltz-2 weights, so its interface score grades
its own family's output. The TNF-alpha campaign measured what that is worth:
Chai-1 against Boltz-2 correlated at Spearman +0.046 on the quantity driving the
submission ranking — no relationship — while reproducing a separate artefact
closely enough to confirm that one was real.

So the question here is not "what does Chai-1 score" but "does the ordering
survive a model that shares no weights with the generator". Three outcomes, and
they imply different submissions:

- ranks agree       -> the ordering is evidence, rank on it
- ranks uncorrelated-> the ordering is a Boltz-2 reading, submit a diverse spread
                       and claim no order
- ranks inverted    -> worse than nothing; the generator's preference is
                       actively misleading here

Binder-to-target only. A global score mixes in the target's own internal
structure, which is the artefact already corrected for on TNF-alpha.
"""

from __future__ import annotations

import csv
import json
import statistics as st
from pathlib import Path
from typing import Annotated

import typer

BINDER_CHAIN_INDEX = 0

app = typer.Typer(add_completion=False)


def spearman(left: list[float], right: list[float]) -> float:
    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        for rank, index in enumerate(order):
            out[index] = float(rank)
        return out

    a, b = ranks(left), ranks(right)
    n = len(a)
    if n < 3:
        return float("nan")
    d2 = sum((x - y) ** 2 for x, y in zip(a, b))
    return 1 - 6 * d2 / (n * (n * n - 1))


def binder_to_target(sample: dict[str, object]) -> float | None:
    """Mean interface score from the binder row to every target chain."""
    matrix = sample.get("per_chain_pair_iptm")
    if not matrix:
        return None
    rows = matrix[0] if isinstance(matrix[0][0], list) else matrix  # type: ignore[index]
    binder_row = rows[BINDER_CHAIN_INDEX]  # type: ignore[index]
    others = [float(v) for i, v in enumerate(binder_row) if i != BINDER_CHAIN_INDEX]
    return st.mean(others) if others else None


def load_chai(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for line in path.open():
        if not line.strip():
            continue
        row = json.loads(line)
        values = [v for s in row.get("samples", []) if (v := binder_to_target(s)) is not None]
        if values:
            out[row["job"]] = st.mean(values)
    return out


@app.command()
def main(
    generator: Annotated[Path, typer.Option(help="pooled BoltzGen metrics CSV")],
    independent: Annotated[Path, typer.Option(help="Chai-1 scores JSONL")],
    out: Annotated[Path, typer.Option(help="per-design comparison CSV")],
) -> None:
    chai = load_chai(independent)
    with generator.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["name"] in chai]
    if len(rows) < 3:
        typer.echo(f"only {len(rows)} designs scored by both; nothing to compare")
        raise typer.Exit(code=1)

    rows.sort(key=lambda r: -float(r["iptm"]))
    boltz = [float(r["iptm"]) for r in rows]
    other = [chai[r["name"]] for r in rows]
    rho = spearman(boltz, other)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "variant", "boltzgen_iptm", "chai_binder_to_target"])
        for r, c in zip(rows, other):
            writer.writerow([r["name"], r["variant"], r["iptm"], f"{c:.3f}"])

    typer.echo(f"{len(rows)} designs scored by both models\n")
    typer.echo(f"{'name':26s} {'variant':11s} {'BoltzGen':>9s} {'Chai-1':>8s}")
    for r, c in zip(rows, other):
        typer.echo(f"{r['name']:26s} {r['variant']:11s} {float(r['iptm']):9.3f} {c:8.3f}")

    typer.echo(f"\nBoltzGen  median {st.median(boltz):.3f}  range {min(boltz):.3f}-{max(boltz):.3f}")
    typer.echo(f"Chai-1    median {st.median(other):.3f}  range {min(other):.3f}-{max(other):.3f}")
    typer.echo(f"\nSpearman, generator vs independent: {rho:+.3f}")
    verdict = (
        "ranks agree - the ordering is evidence" if rho >= 0.5
        else "ranks inverted - the generator's preference is misleading" if rho <= -0.3
        else "no relationship - the ordering is a generator reading, not a measurement"
    )
    typer.echo(f"verdict: {verdict}")
    typer.echo(f"\nwrote {out}")


if __name__ == "__main__":
    app()
