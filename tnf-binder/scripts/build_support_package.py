"""Assemble the supporting-information package that goes up with the submission.

The competition FAQ invites metrics, structure models and methodology alongside
the CSV, and states that for Tracks 2 and 3 the submitted information is pooled
with the designs when selections are made. So the package has to stand on its
own as a technical record: every number carries its denominator and its sampling
depth, and nothing is a point estimate where the measurement was a distribution.

It is assembled rather than written by hand so that it can be regenerated when a
measurement changes -- three of the twenty rows are waiting on a deeper mouse
fold at the time of writing.
"""

from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path
from typing import Annotated

import typer

# Columns carried over from the v8 evidence file, which holds the structural and
# novelty measurements. Renamed where the original name hid the denominator.
CARRY = {
    "interface_area_a2": "interface_area_a2_model0",
    "clashes": "clashes_model0",
    "cysteines": "cysteines",
    "fold_similarity": "fold_similarity_class",
    "proteinbase_novelty": "proteinbase_novelty_prior_round",
    "delta_q_floor": "delta_charge_at_ph6",
    "his_count": "his_total",
    "ideal_window_his": "window_his_3_samples",
}

app = typer.Typer(add_completion=False)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_metrics(v9_evidence: Path, v8_evidence: Path, deep_geometry: list[Path]) -> list[dict[str, object]]:
    v8 = {r["name"]: r for r in read_csv(v8_evidence)}
    geometry = aggregate_geometry(deep_geometry)

    rows: list[dict[str, object]] = []
    for row in read_csv(v9_evidence):
        name = row["name"]
        source = v8.get(name, {})
        record: dict[str, object] = {
            "submission_rank": row["rank"],
            "name": name,
            "design_id": row["design"],
            "chassis": "helix_0" if row["design"].startswith("helix0") else "free_1",
            "length_aa": row["length"],
            "molecule_class": "protein",
            "switch_his_n": row["switch_his"],
            "switch_his_positions": row["switch_positions"],
            "switch_measured_at_20_samples": row["measured_at_20_samples"],
            "his_in_window_rates_n20": row["all_his_rates_n20"],
            "mouse_retention": row["mouse_retention"],
            "mouse_retention_samples": row["retention_samples"],
            "human_binder_to_target_iptm": row["human_iptm"],
        }
        for original, renamed in CARRY.items():
            record[renamed] = source.get(original, "")
        record |= geometry.get(row["design"], {})
        rows.append(record)
    return rows


def aggregate_geometry(paths: list[Path]) -> dict[str, dict[str, object]]:
    """Per-design interface geometry over the twenty-sample fold.

    The campaign's own single-model numbers are kept alongside these because one
    submitted row was nearly discarded on a model_0 that read 232 A2 against a
    twenty-sample median of 825.
    """
    # Rows are per diffusion sample, named {CAMPAIGN}__{design}_model_{n}. Both
    # the campaign prefix and the sample suffix have to come off, or every sample
    # is aggregated as its own design and every statistic is a single model again.
    by_design: dict[str, list[dict[str, str]]] = {}
    for path in paths:
        if not path.exists():
            continue
        for row in read_csv(path):
            label = row.get("design") or row.get("name") or ""
            design = re.sub(r"_model_\d+$", "", re.sub(r"^[A-Za-z0-9]+__", "", label))
            by_design.setdefault(design, []).append(row)

    out: dict[str, dict[str, object]] = {}
    for design, samples in by_design.items():
        areas = sorted(float(s["interface_area"]) for s in samples if s.get("interface_area"))
        clashfree = sum(1 for s in samples if s.get("clashes") in {"0", "0.0"})
        if not areas:
            continue
        out[design] = {
            "interface_area_median_n20": round(areas[len(areas) // 2], 1),
            "interface_area_min_n20": round(areas[0], 1),
            "clash_free_samples_n20": f"{clashfree}/{len(samples)}",
        }
    return out


@app.command()
def main(
    v9_evidence: Annotated[Path, typer.Option(help="v9 evidence CSV")],
    v8_evidence: Annotated[Path, typer.Option(help="v8 evidence CSV, for structural metrics")],
    deep_geometry: Annotated[list[Path], typer.Option(help="per-sample geometry CSVs at n=20")],
    structures: Annotated[Path, typer.Option(help="directory of per-design CIF/PDB")],
    docs: Annotated[list[Path], typer.Option(help="methodology documents to include")],
    out_dir: Annotated[Path, typer.Option(help="package directory to build")],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = build_metrics(v9_evidence, v8_evidence, deep_geometry)

    metrics = out_dir / "metrics.csv"
    with metrics.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    doc_dir = out_dir / "methodology"
    doc_dir.mkdir(exist_ok=True)
    copied = []
    for doc in docs:
        if doc.exists():
            shutil.copy2(doc, doc_dir / doc.name)
            copied.append(doc.name)

    struct_dir = out_dir / "structures"
    struct_dir.mkdir(exist_ok=True)
    n_struct = 0
    manifest: list[dict[str, str]] = []
    for row in rows:
        name = str(row["name"])
        for suffix in (".cif", ".pdb"):
            source = structures / f"{name}{suffix}"
            if source.exists():
                shutil.copy2(source, struct_dir / f"rank{int(row['submission_rank']):02d}_{name}{suffix}")
                n_struct += 1
        manifest.append({
            "submission_rank": str(row["submission_rank"]),
            "name": name,
            "design_id": str(row["design_id"]),
            "file": f"rank{int(row['submission_rank']):02d}_{name}.cif",
        })
    with (struct_dir / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    typer.echo(f"metrics.csv      {len(rows)} designs x {len(rows[0])} columns")
    typer.echo(f"methodology/     {len(copied)} documents: {', '.join(copied)}")
    typer.echo(f"structures/      {n_struct} files + manifest")
    typer.echo(f"built {out_dir}")


if __name__ == "__main__":
    app()
