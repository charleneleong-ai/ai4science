"""Measure what the available in-silico scores are actually worth.

The organisers' own post-mortem says their ranking metrics did not predict
binding, and an independent re-analysis reports that confidence calibration
inverts above 0.6. Both claims are checkable: the Proteinbase corpus publishes
predicted scores and wet-lab outcomes for the same designs. This recomputes
them rather than taking either on trust, because the result decides how 20
submission slots get ranked.

Pooled AUROC across targets is the wrong headline -- easy targets inflate it --
so within-target AUROC is reported alongside, weighted by the number of
comparisons each target supports.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Annotated

import typer

SCORES = (
    "boltz2_iptm",
    "boltz2_pdockq",
    "boltz2_pdockq2",
    "boltz2_lis",
    "boltz2_plddt",
    "boltz2_complex_plddt",
    "boltz2_complex_iplddt",
    "boltz2_complex_pde",
    "esmfold_plddt",
    "seqidentity",
    "proteinmpnn_score",
)

# Adaptyv's post-mortem reports expression is predicted by composition where
# structural scores fail, so these are computed from the sequence as controls.
COMPOSITION = ("glu_fraction", "lys_fraction")

app = typer.Typer(add_completion=False)


def auroc(positives: list[float], negatives: list[float]) -> tuple[float, int]:
    """Mann-Whitney AUROC with ties counted as half, plus the comparison count."""
    if not positives or not negatives:
        return float("nan"), 0
    wins = 0.0
    for p in positives:
        for n in negatives:
            wins += 1.0 if p > n else 0.5 if p == n else 0.0
    return wins / (len(positives) * len(negatives)), len(positives) * len(negatives)


def load_designs(path: Path) -> list[dict[str, object]]:
    designs: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            raw = row.get("evaluations")
            if not raw:
                continue
            record: dict[str, object] = {"name": row["name"], "method": row.get("designMethod") or "unknown"}
            sequence = row.get("sequence") or ""
            if sequence:
                record["glu_fraction"] = sequence.count("E") / len(sequence)
                record["lys_fraction"] = sequence.count("K") / len(sequence)
            for evaluation in json.loads(raw):
                record[evaluation["metric"]] = evaluation["value"]
                if evaluation["type"] == "experimental" and evaluation.get("target"):
                    record["target"] = evaluation["target"]
            designs.append(record)
    return designs


def as_float(value: object) -> float | None:
    # seqidentity is published as a nested object rather than a number.
    if isinstance(value, dict):
        value = value.get("value")
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def label(record: dict[str, object], metric: str) -> bool | None:
    value = record.get(metric)
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in {"true", "false"}:
        return value.lower() == "true"
    return None


def constant_fields(designs: list[dict[str, object]], names: tuple[str, ...]) -> dict[str, float]:
    """Fields with a single distinct value are placeholders, not weak predictors.

    Reporting AUROC 0.500 for one of these would say "no signal" when the truth
    is "no data", so they are named and excluded.
    """
    frozen: dict[str, float] = {}
    for name in names:
        values = {v for d in designs if (v := as_float(d.get(name))) is not None}
        if len(values) == 1:
            frozen[name] = values.pop()
    return frozen


def reliability(pairs: list[tuple[float, bool]], edges: tuple[float, ...]) -> list[tuple[str, int, float]]:
    """Observed hit rate per predicted-score band, to test the inversion claim."""
    out: list[tuple[str, int, float]] = []
    for lo, hi in zip(edges, edges[1:]):
        band = [hit for score, hit in pairs if lo <= score < hi]
        if band:
            out.append((f"{lo:.1f}-{hi:.1f}", len(band), sum(band) / len(band)))
    return out


@app.command()
def main(
    corpus: Annotated[Path, typer.Option(help="proteinbase_latest.csv")],
    out: Annotated[Path, typer.Option(help="summary CSV")],
    outcome: Annotated[str, typer.Option(help="experimental metric to predict: binding or expressed")] = "binding",
) -> None:
    designs = load_designs(corpus)
    labelled = [d for d in designs if label(d, outcome) is not None]
    positives = sum(1 for d in labelled if label(d, outcome))
    typer.echo(
        f"{len(designs)} designs, {len(labelled)} with a '{outcome}' label,"
        f" {positives} positive ({positives / max(len(labelled), 1):.1%})"
    )

    frozen = constant_fields(designs, SCORES + COMPOSITION)
    if frozen:
        typer.echo("\nconstant fields -- no data, excluded from the table:")
        for name, value in frozen.items():
            typer.echo(f"  {name:24s} every row = {value}")

    by_target: dict[str, list[dict[str, object]]] = {}
    for record in labelled:
        by_target.setdefault(str(record.get("target", "unknown")), []).append(record)

    typer.echo(f"\nper-target {outcome} rates (n >= 20):")
    for target, records in sorted(by_target.items(), key=lambda kv: -len(kv[1])):
        if len(records) < 20:
            continue
        hits = sum(1 for r in records if label(r, outcome))
        typer.echo(f"  {target:22s} n={len(records):4d}  binders={hits:4d}  {hits / len(records):6.1%}")

    rows: list[dict[str, object]] = []
    typer.echo(f"\n{'score':24s} {'pooled':>8s} {'n_pool':>8s} {'within':>8s} {'targets':>8s}")
    for score in [s for s in SCORES + COMPOSITION if s not in frozen]:
        positives = [v for r in labelled if label(r, outcome) and (v := as_float(r.get(score))) is not None]
        negatives = [v for r in labelled if not label(r, outcome) and (v := as_float(r.get(score))) is not None]
        pooled, n_pool = auroc(positives, negatives)

        weighted, weight, used = 0.0, 0, 0
        for records in by_target.values():
            pos = [v for r in records if label(r, outcome) and (v := as_float(r.get(score))) is not None]
            neg = [v for r in records if not label(r, outcome) and (v := as_float(r.get(score))) is not None]
            value, count = auroc(pos, neg)
            if count >= 20:
                weighted += value * count
                weight += count
                used += 1
        within = weighted / weight if weight else float("nan")

        rows.append(
            {"score": score, "pooled_auroc": round(pooled, 3), "n_pooled": n_pool,
             "within_target_auroc": round(within, 3), "targets": used}
        )
        typer.echo(f"{score:24s} {pooled:8.3f} {n_pool:8d} {within:8.3f} {used:8d}")

    pairs = [
        (v, bool(label(r, outcome)))
        for r in labelled
        if (v := as_float(r.get("boltz2_iptm"))) is not None
    ]
    typer.echo(f"\nreliability of boltz2_iptm against {outcome} (does the curve invert above 0.6?):")
    for band, n, rate in reliability(pairs, (0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.01)):
        typer.echo(f"  iptm {band:10s} n={n:4d}  observed binding {rate:6.1%}")

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    typer.echo(f"\nwrote {out}")


if __name__ == "__main__":
    app()
