"""Re-rank the submission on agreement between methods rather than on one method.

Three independent measurements now exist per design-position: PROPKA's in-window
rate over twenty samples, pKAI's in-window rate over the same twenty, and the
bound-versus-free pKa linkage that says whether the titration is coupled to
binding at all.

They are combined as tiers rather than as a single hard gate, so that a check
which turns out to reject everything degrades the ranking gracefully instead of
flattening it. A design is ordered by how many histidines reach the strongest
tier it can, then by the next tier down, and only then by the original
objectives.

    tier A  both pKa methods in-window on >=90% of samples AND coupled to binding
    tier B  both pKa methods agree, coupling absent or too weak to call
    tier C  one method only

Tier C is reported but does not promote: a switch only one method believes in is
exactly the claim this exercise was built to stop making.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import typer

SWITCH_RATE = 0.90
# Matches COUPLING_FLOOR in pka_binding_coupling.py: pKa(complex) - pKa(free)
# must be this far below zero before protonation is called binding-weakening.
COUPLING_FLOOR = -0.3

app = typer.Typer(add_completion=False)


def strip_campaign(design: str) -> str:
    return re.sub(r"^[A-Za-z0-9]+__", "", design)


def read_rates(path: Path, column: str) -> dict[tuple[str, int], float]:
    out: dict[tuple[str, int], float] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            key = (strip_campaign(row["design"]), int(row["position"]))
            if column == "in_window_samples":
                hit, total = row[column].split("/")
                out[key] = int(hit) / int(total)
            else:
                out[key] = float(row[column])
    return out


@dataclass
class Scored:
    name: str
    design: str
    sequence: str
    retention: float
    retention_n: int
    iptm: float
    tier_a: list[int] = field(default_factory=list)
    tier_b: list[int] = field(default_factory=list)
    tier_c: list[int] = field(default_factory=list)

    @property
    def rank_key(self) -> tuple[int, int, float, float]:
        return (-len(self.tier_a), -len(self.tier_b), -self.retention, -self.iptm)


@app.command()
def main(
    evidence: Annotated[Path, typer.Option(help="current v9 evidence CSV")],
    submitted: Annotated[Path, typer.Option(help="current v9 submission CSV")],
    propka: Annotated[Path, typer.Option(help="pka_robustness_deep.csv")],
    pkai: Annotated[Path, typer.Option(help="pka_pkai_summary.csv")],
    coupling: Annotated[Path, typer.Option(help="pka_coupling_summary.csv")],
    out: Annotated[Path, typer.Option(help="re-ranked submission CSV")],
    out_evidence: Annotated[Path, typer.Option(help="re-ranked evidence CSV")],
    molecule_class: Annotated[str, typer.Option(help="class value the uploader accepts")] = "single_chain",
) -> None:
    propka_rate = read_rates(propka, "in_window_samples")
    pkai_rate = read_rates(pkai, "in_window_rate")
    delta = read_rates(coupling, "delta_mean")

    with submitted.open() as handle:
        sequences = {r["name"]: r["sequence"] for r in csv.DictReader(handle)}

    designs: list[Scored] = []
    with evidence.open() as handle:
        for row in csv.DictReader(handle):
            design = row["design"]
            scored = Scored(
                name=row["name"],
                design=design,
                sequence=sequences[row["name"]],
                retention=float(row["mouse_retention"] or 0),
                retention_n=int(row["retention_samples"] or 0),
                iptm=float(row["human_iptm"] or 0),
            )
            positions = {p for (d, p) in propka_rate if d == design} | {
                p for (d, p) in pkai_rate if d == design
            }
            for position in sorted(positions):
                key = (design, position)
                by_propka = propka_rate.get(key, 0.0) >= SWITCH_RATE
                by_pkai = pkai_rate.get(key, 0.0) >= SWITCH_RATE
                coupled = delta.get(key)
                if by_propka and by_pkai and coupled is not None and coupled <= COUPLING_FLOOR:
                    scored.tier_a.append(position)
                elif by_propka and by_pkai:
                    scored.tier_b.append(position)
                elif by_propka or by_pkai:
                    scored.tier_c.append(position)
            designs.append(scored)

    designs.sort(key=lambda d: d.rank_key)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "sequence", "molecule_class"])
        for d in designs:
            writer.writerow([d.name, d.sequence, molecule_class])

    with out_evidence.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([
            "rank", "name", "design", "tier_a_both_methods_and_coupled", "tier_b_both_methods",
            "tier_c_one_method_only", "mouse_retention", "retention_samples", "human_iptm",
        ])
        for i, d in enumerate(designs, 1):
            writer.writerow([
                i, d.name, d.design,
                " ".join(f"H{p}" for p in d.tier_a) or "none",
                " ".join(f"H{p}" for p in d.tier_b) or "none",
                " ".join(f"H{p}" for p in d.tier_c) or "none",
                f"{d.retention:.3f}", d.retention_n, f"{d.iptm:.3f}",
            ])

    typer.echo(f"{'#':>3s} {'name':9s} {'design':24s} {'A':>2s} {'B':>2s} {'C':>2s} {'retention':>10s}")
    for i, d in enumerate(designs, 1):
        typer.echo(
            f"{i:3d} {d.name:9s} {d.design:24s} {len(d.tier_a):2d} {len(d.tier_b):2d}"
            f" {len(d.tier_c):2d} {d.retention:10.3f}"
        )
    typer.echo(f"\nwrote {out} and {out_evidence}")


if __name__ == "__main__":
    app()
