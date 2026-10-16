"""Rebuild the Challenge 2 submission against the ranked objectives.

Three defects in v8, found while reading the challenge FAQ for Challenge 3:

1. `molecule_class` was `single_chain`. The permitted values are "protein",
   "nanobody", "scfv", "fab_kappa" and "fab_lambda", and the FAQ states that
   length and class annotation are *hard* filters applied to every design.
2. The rows were emitted in name order, not ranked order. The FAQ asks for a CSV
   "ordered by how you would rank your molecules (top row higher)", and for
   Tracks 2 and 3 the submitted ordering and information are what the selection
   workflow reads.
3. The ranking that did exist counted any histidine landing in the pH window on
   three samples. At twenty samples His21, His35 and His38 have in-window rates
   of 0.56, 0.43 and 0.56 pooled, so they are not switches; only His67 (0.94),
   His72 (0.84) and His70 (0.83) hold up, and the rate is design-specific rather
   than positional.

The challenge ranks pH-selectivity first, mouse cross-reactivity second and
affinity third, so that is the sort order here, each term measured at the
deepest sampling available and labelled with its denominator.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer

# The two official sources disagree. FAQ 5 and the challenge page both enumerate
# the permitted values as "protein", "nanobody", "scfv", "fab_kappa",
# "fab_lambda"; the downloadable submission template labels its single-chain
# example "single_chain". Three of the template's four values match the FAQ
# exactly and only the single-chain one differs, which reads like a rename the
# template did not pick up -- but that is an inference, so both variants are
# written and the choice is made at upload.
MOLECULE_CLASS = "protein"
MOLECULE_CLASS_TEMPLATE = "single_chain"

# The template is LF-terminated and unquoted; csv.writer defaults to CRLF.
LINE_TERMINATOR = "\n"

# A histidine counts as a switch if its predicted pKa sits in the 6.0-7.6 window
# on at least this fraction of twenty samples. Declared post hoc on the full
# sensitivity curve (1.00 selects 0 designs, 0.95 selects 2, 0.90 selects 5,
# 0.85 selects 7, 0.80 selects 9, 0.75 selects 10) rather than as a single cut.
SWITCH_RATE = 0.90

app = typer.Typer(add_completion=False)


# Designs whose twenty-sample mouse fold was measured outside the main deep run:
# the free_1 pair in its own deep run, and the three top-tier designs folded
# 2026-10-17 specifically because the submission order turned on their retention
# and they only had three samples. All three three-sample values were optimistic
# (0.874 -> 0.802, 0.824 -> 0.735, 0.853 -> 0.860), which reordered ranks 3-5.
EXTRA_RETENTION = {
    "free1_T73H": (0.466, 20),
    "free1_parent": (0.482, 20),
    "helix0_E67H_E70H": (0.860, 20),
    "helix0_E67H_R72H": (0.802, 20),
    "helix0_R21H_E67H_E70H": (0.735, 20),
}
EXTRA_IPTM = {"free1_T73H": 0.434, "free1_parent": 0.424}


@dataclass
class Design:
    name: str
    sequence: str
    design: str
    switches: tuple[int, ...]
    switch_rates: dict[int, float]
    provisional_switches: int
    measured_at_20: bool
    retention: float | None
    retention_n: int
    iptm: float

    @property
    def switch_count(self) -> int:
        return len(self.switches) if self.measured_at_20 else self.provisional_switches

    @property
    def rank_key(self) -> tuple[int, int, float, float]:
        # Objective order: pH-selectivity, then mouse cross-reactivity, then
        # affinity. Within a switch count, a design measured at twenty samples
        # outranks one whose count is still provisional on three -- an unmeasured
        # design must not be scored as if it had failed, and must not outrank a
        # design that actually passed.
        return (
            -self.switch_count,
            0 if self.measured_at_20 else 1,
            -(self.retention or 0.0),
            -self.iptm,
        )


def load_switch_rates(paths: list[Path]) -> dict[str, dict[int, float]]:
    rates: dict[str, dict[int, float]] = {}
    for path in paths:
        with path.open() as handle:
            for row in csv.DictReader(handle):
                hit, total = row["in_window_samples"].split("/")
                # Campaign prefixes (A__, F__) tag which run produced the fold;
                # the design itself is the same molecule.
                design = re.sub(r"^[A-Z]__", "", row["design"])
                rates.setdefault(design, {})[int(row["position"])] = int(hit) / int(total)
    return rates


def load_retention(path: Path) -> dict[str, tuple[float, int]]:
    out: dict[str, tuple[float, int]] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 6 or not parts[0].startswith("helix0"):
            continue
        out[parts[0]] = (float(parts[5]), int(parts[2]))
    return out


def load_evidence(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


@app.command()
def main(
    evidence: Annotated[Path, typer.Option(help="v8 evidence CSV")],
    submitted: Annotated[Path, typer.Option(help="v8 submission CSV, for the sequences")],
    deep_pka: Annotated[list[Path], typer.Option(help="pKa robustness CSVs at n=20")],
    retention: Annotated[Path, typer.Option(help="20-sample retention table")],
    out: Annotated[Path, typer.Option(help="v9 submission CSV")],
    out_evidence: Annotated[Path, typer.Option(help="v9 evidence CSV")],
) -> None:
    rates = load_switch_rates(deep_pka)
    retained = load_retention(retention) | EXTRA_RETENTION
    sequences = {r["name"]: r["sequence"] for r in load_evidence(submitted)}

    designs: list[Design] = []
    for row in load_evidence(evidence):
        name = row["name"]
        design = row["design"]
        per_position = rates.get(design, {})
        switches = tuple(sorted(p for p, rate in per_position.items() if rate >= SWITCH_RATE))
        # A design with no histidine at all is measured by construction.
        measured = bool(per_position) or int(row["his_count"] or 0) == 0
        value = retained.get(design)
        three_sample = row.get("mouse_retained") or row.get("cross_reactivity") or ""
        if value is not None:
            retention_value, retention_n = value
        else:
            try:
                retention_value, retention_n = float(three_sample), 3
            except ValueError:
                retention_value, retention_n = None, 0
        designs.append(
            Design(
                name=name,
                sequence=sequences[name],
                design=design,
                switches=switches,
                switch_rates=per_position,
                provisional_switches=int(row["ideal_window_his"] or 0),
                measured_at_20=measured,
                retention=retention_value,
                retention_n=retention_n,
                iptm=EXTRA_IPTM.get(design, float(row["binder_to_target_iptm"] or 0)),
            )
        )

    designs.sort(key=lambda d: d.rank_key)

    out.parent.mkdir(parents=True, exist_ok=True)
    for path, molecule_class in (
        (out, MOLECULE_CLASS),
        (out.with_name(f"{out.stem}_{MOLECULE_CLASS_TEMPLATE}{out.suffix}"), MOLECULE_CLASS_TEMPLATE),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator=LINE_TERMINATOR)
            writer.writerow(["name", "sequence", "molecule_class"])
            for d in designs:
                writer.writerow([d.name, d.sequence, molecule_class])

    with out_evidence.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator=LINE_TERMINATOR)
        writer.writerow(
            ["rank", "name", "design", "length", "switch_his", "measured_at_20_samples",
             "switch_positions", "all_his_rates_n20", "mouse_retention", "retention_samples",
             "human_iptm"]
        )
        for i, d in enumerate(designs, 1):
            writer.writerow([
                i, d.name, d.design, len(d.sequence), d.switch_count,
                "yes" if d.measured_at_20 else "no, provisional on 3 samples",
                " ".join(f"H{p}" for p in d.switches) or "none",
                " ".join(f"H{p}={r:.2f}" for p, r in sorted(d.switch_rates.items())) or "not folded at n=20",
                f"{d.retention:.3f}" if d.retention is not None else "",
                d.retention_n, f"{d.iptm:.3f}",
            ])

    typer.echo(f"{'#':>3s} {'name':9s} {'design':24s} {'sw':>3s} {'retention':>12s} {'ipTM':>6s}")
    for i, d in enumerate(designs, 1):
        flag = "" if d.retention_n == 20 else f"  ret n={d.retention_n}"
        if not d.measured_at_20:
            flag += "  SWITCH PROVISIONAL"
        ret = f"{d.retention:.3f}" if d.retention is not None else "-"
        typer.echo(f"{i:3d} {d.name:9s} {d.design:24s} {d.switch_count:3d} {ret:>12s} {d.iptm:6.3f}{flag}")
    typer.echo(f"\nwrote {out} and {out_evidence}")


if __name__ == "__main__":
    app()
