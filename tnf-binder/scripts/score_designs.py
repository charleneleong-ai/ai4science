"""Score a batch of predicted binder/TNF-alpha complexes on the objectives the generator skips.

`tnf-binder verify` is the one-structure version of this and shares every threshold through
`InterfaceVerifier`; this is the batch form, which is what a campaign actually needs: one CSV
row per design, ranked, with a bad file costing one row rather than the run.

Everything is derived from per-residue buried surface area (unbound SASA minus bound SASA),
not from atom- or residue-pair counts. Pair counting was tried first and was badly misleading:
it scaled with interpenetration, so clashing poly-alanine junk out-scored real designs, and the
per-residue tallies were not even physically plausible. Buried area is a direct physical
quantity and degrades gracefully.

interface_area    target-side buried area, A^2. A real minibinder interface is ~600-1000.
epitope_focus     share of the target's buried area contributed by the intended epitope.
cross_reactivity  1 - share of buried area on human/mouse-divergent positions. Objective 2.

There is deliberately no geometric pH-switch column. The earlier one (His burial x standoff
from a conserved cation) was refuted by PROPKA on real designs -- burial and cation proximity
both correlate with the WRONG sign, see `verifiers/pka_switch.py`. Objective 3 is measured by
`his_pka_scan.py` instead.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.assess import assess
from tnf_binder.tables import write_rows
from tnf_binder.target import LOBE, POCKET, TnfTarget

app = typer.Typer(add_completion=False)
# stdout is the CSV itself, so every message goes to stderr -- `script dir > scored.csv` is
# how this is used, and a wrapped rich table in the middle of it would corrupt the file.
console = Console(stderr=True)

EPITOPES = {"pocket": POCKET, "lobe": LOBE}
COLUMNS = ["design", "binder_len", "interface_area", "epitope_focus", "cross_reactivity",
           "divergent_frac", "clashes", "flag", "sequence"]


def score(path: Path, epitope: tuple[int, ...], divergent: frozenset[int]) -> dict[str, object]:
    """One CSV row per complex. Every verifier's reason lands in `flag`, not just the
    interface one: the top of a 40-design shortlist was poly-alanine with a perfectly good
    interface, and only composition saw it."""
    assessment = assess(path, epitope=epitope, divergent=divergent)
    if not assessment.is_resolved:
        return {"design": path.stem, "flag": "chain_classification_failed"}
    design, metrics = assessment.design, assessment.metrics
    return {
        "design": path.stem,
        "binder_len": design.length,
        "interface_area": round(metrics.area, 1),
        "epitope_focus": round(metrics.epitope_focus, 4),
        "cross_reactivity": round(metrics.cross_reactivity, 4),
        "divergent_frac": round(metrics.divergent_fraction, 4),
        "clashes": metrics.clashes,
        "flag": ",".join(assessment.failures),
        "sequence": design.sequence,
    }


def prefer_refold(paths: list[Path]) -> list[Path]:
    """One path per design stem, preferring the refold over the design-step pose.

    A pipeline output dir holds BOTH for the same stem and they are not interchangeable: the
    design step places grafted scaffolds only approximately (a nanobody measured 524 atom pairs
    under 2.2 A, interpenetrating at backbone level), while the refold is the model's own
    physically sensible re-prediction. Scoring a parent dir silently produced two contradictory
    rows per design.
    """
    by_stem: dict[str, Path] = {}
    for path in paths:
        previous = by_stem.get(path.stem)
        if previous is None or ("refold" in path.parts and "refold" not in previous.parts):
            by_stem[path.stem] = path
    return sorted(by_stem.values())


@app.command()
def main(
    structures: Path = typer.Argument(..., help="directory of predicted complexes"),
    out: Path | None = typer.Option(None, "--out", help="also write the table here"),
    epitope: str = typer.Option("pocket", "--epitope",
                                help=f"which epitope focus is measured against: "
                                     f"{', '.join(EPITOPES)}"),
) -> None:
    """Score every complex under a directory, best interface first."""
    if epitope not in EPITOPES:
        typer.echo(f"unknown epitope {epitope!r}; choose from {', '.join(EPITOPES)}", err=True)
        raise typer.Exit(code=2)
    found = sorted(p for p in structures.rglob("*") if p.suffix in (".cif", ".pdb", ".mmcif"))
    if not found:
        typer.echo(f"no structures under {structures}", err=True)
        raise typer.Exit(code=2)

    paths = prefer_refold(found)
    if len(paths) != len(found):
        console.print(f"[dim]{len(found)} files collapsed to {len(paths)} designs "
                      f"(preferred refold where both exist)[/]")

    divergent = TnfTarget().divergent
    rows = []
    for path in paths:
        try:
            rows.append(score(path, EPITOPES[epitope], divergent))
        except Exception as exc:  # one bad file must not sink the batch
            rows.append({"design": path.stem, "flag": f"scoring_error:{type(exc).__name__}"})
    rows.sort(key=lambda r: (bool(r.get("flag")), -float(r.get("epitope_focus") or 0),
                             -float(r.get("interface_area") or 0)))

    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    if out:
        write_rows(out, rows, COLUMNS)
        console.print(f"\nwrote [cyan]{out}[/]")


if __name__ == "__main__":
    app()
