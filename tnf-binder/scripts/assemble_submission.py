"""Assemble ranked per-track shortlists from every measurement, applying the hard gates.

`make_final_csv.py` merges these into the one file that is uploaded; this stage is where the
gates run and where the metrics stay visible beside each sequence.

Ranking follows the competition's own order -- pH-selectivity, then mouse cross-reactivity,
then affinity -- and its statement that a weak-but-clearly-pH-selective binder is worth more
than a high-affinity one that is not selective.

Interface confidence is written as a MODEL-DEPENDENT value, never as truth: Boltz-2 and
BoltzGen's folding step disagreed by 2-3x on identical structures, and AF2 scored its own
optimised designs at 0.92 where an independent model scored them 0.21.

The gates below drop artefacts rather than down-ranking them. That distinction is load-bearing:
the top-scoring construct of one campaign was verbatim ubiquitin, and the top of a 40-design
shortlist was poly-alanine. Ranking cannot rescue either.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.design import TnfDesign
from tnf_binder.submission import MOLECULE_CLASS, row_problem, write_submission
from tnf_binder.tables import number, read_rows
from tnf_binder.verifiers import CompositionVerifier

app = typer.Typer(add_completion=False)
console = Console()

# Shortlist sizes per track. The upload itself is capped at 20 across all tracks, which
# make_final_csv.py applies; these only bound how much each ranking is allowed to carry.
LIMITS = {1: 40, 2: 20, 3: 20}

COLUMNS = ["name", "sequence", "molecule_class", "boltz_iptm_best", "boltz_iptm_mean",
           "boltz_ipde", "n_his_ideal", "ph_switch_score", "pka_detail",
           "cross_reactivity", "epitope_focus", "notes"]

NOTE = ("epitope=conserved_TNFR2_site; pH switch = PROPKA-verified His pKa in 6.0-7.6; "
        "iptm is Boltz-2 binder-to-target (BoltzGen folding gave ~0.2 on the same designs "
        "-- models disagree, both reported)")

# Objective order, from the competition's own ranking.
TRACK_ORDER: dict[int, Callable[[dict[str, object]], tuple[float, ...]]] = {
    1: lambda r: (-r["boltz_iptm_mean"], -r["ph_switch_score"]),
    2: lambda r: (-(r["cross_reactivity"] or 0), -r["boltz_iptm_mean"]),
    3: lambda r: (-r["ph_switch_score"], -r["boltz_iptm_mean"]),
}


def gate_reason(name: str, sequence: str, novel: bool | None,
                geometry: dict[str, str]) -> str | None:
    """Why a construct is disqualified, or None if it survives every hard gate."""
    if problem := row_problem({"name": name, "sequence": sequence,
                               "molecule_class": MOLECULE_CLASS}):
        return problem
    if novel is False:
        return "MEMORISED (novelty gate)"
    if novel is None:
        return "no novelty verdict"
    composition = CompositionVerifier().verify(TnfDesign(name, sequence))
    if not composition.trust:
        return composition.reason
    flag = (geometry.get("flag") or "").strip()
    return f"structural flag [{flag}]" if flag else None


def ph_switch_terms(row: dict[str, str], pka_by_position: dict[tuple[str, str], dict[str, str]],
                    ) -> tuple[int, int, list[str]]:
    """PROPKA-verified in-window histidines for one construct."""
    ideal = usable = 0
    detail = []
    for token in (row.get("positions") or "").split():
        record = pka_by_position.get((row["design"], token))
        if not record:
            continue
        detail.append(f"{token}:pKa={record['pka']}")
        ideal += record["verdict"] == "ideal"
        usable += record["verdict"] == "usable"
    return ideal, usable, detail


def novelty_lookup(novelty_rows: list[dict[str, str]]) -> Callable[[str], bool | None]:
    """spec -> bool, tolerant of the run/chain suffixes the various tools append."""
    if not novelty_rows:
        return lambda _spec: True
    verdict = {r["design"]: r["novel"].strip().lower() == "true" for r in novelty_rows}

    def check(spec: str) -> bool | None:
        if spec in verdict:
            return verdict[spec]
        return next((v for k, v in verdict.items() if spec in k or k in spec), None)

    return check


def build_record(row: dict[str, str], pka_by_position: dict[tuple[str, str], dict[str, str]],
                 geometry: dict[str, str]) -> dict[str, object]:
    ideal, usable, detail = ph_switch_terms(row, pka_by_position)
    return {
        "name": row["spec"],
        "sequence": (row.get("sequence") or "").strip().upper(),
        "boltz_iptm_best": round(number(row, "best"), 3),
        "boltz_iptm_mean": round(number(row, "mean"), 3),
        "boltz_ipde": round(number(row, "ipde"), 2),
        "n_his_ideal": ideal,
        "ph_switch_score": round(ideal + 0.4 * usable, 2),
        "pka_detail": ";".join(detail),
        "cross_reactivity": round(number(geometry, "cross_reactivity"), 4) if geometry else "",
        "epitope_focus": round(number(geometry, "epitope_focus"), 4) if geometry else "",
    }


def dedupe(records: list[dict[str, object]]) -> list[dict[str, object]]:
    seen, unique = set(), []
    for record in sorted(records, key=lambda r: -r["boltz_iptm_mean"]):
        if record["sequence"] not in seen:
            seen.add(record["sequence"])
            unique.append(record)
    return unique


def write_track(path: Path, records: list[dict[str, object]], track: int, prefix: str) -> None:
    write_submission(path, COLUMNS, [
        [f"{prefix}_t{track}_{i:02d}", r["sequence"], MOLECULE_CLASS, r["boltz_iptm_best"],
         r["boltz_iptm_mean"], r["boltz_ipde"], r["n_his_ideal"], r["ph_switch_score"],
         r["pka_detail"], r["cross_reactivity"], r["epitope_focus"],
         f"src={r['name']}; {NOTE}"]
        for i, r in enumerate(records, 1)])


@app.command()
def main(
    campaign: Path = typer.Option(..., "--campaign", help="measured rows from analyze_campaign"),
    outdir: Path = typer.Option(..., "--outdir"),
    novelty: Path | None = typer.Option(None, "--novelty"),
    geometry: Path | None = typer.Option(None, "--geometry", help="table from score_designs"),
    pka: list[Path] = typer.Option([], "--pka", help="pKa scan CSVs, repeatable"),
    prefix: str = typer.Option("tnf", "--prefix"),
) -> None:
    """Gate every construct, then write one ranked shortlist per objective."""
    rows = read_rows(campaign)
    pka_by_position = {(r["parent"], f"{r['chain']}{r['resnum']}"): r
                       for path in pka for r in read_rows(path)}
    geometry_by_design = {r["design"]: r for r in read_rows(geometry)}
    is_novel = novelty_lookup(read_rows(novelty))
    console.print(f"{len(rows)} constructs from campaign")

    kept, dropped = [], []
    for row in rows:
        sequence = (row.get("sequence") or "").strip().upper()
        geometry_row = geometry_by_design.get(row["spec"], {})
        if reason := gate_reason(row["spec"], sequence, is_novel(row["spec"]), geometry_row):
            dropped.append(f"{row['spec']}: {reason}")
        else:
            kept.append(build_record(row, pka_by_position, geometry_row))

    console.print(f"gates: kept [bold]{len(kept)}[/], dropped [bold]{len(dropped)}[/]")
    for message in dropped[:12]:
        console.print(f"  [dim]{message}[/]")

    unique = dedupe(kept)
    console.print(f"unique sequences: {len(unique)}")
    outdir.mkdir(parents=True, exist_ok=True)
    for track, limit in LIMITS.items():
        selected = sorted(unique, key=TRACK_ORDER[track])[:limit]
        path = outdir / f"{prefix}_track{track}.csv"
        write_track(path, selected, track, prefix)
        console.print(f"track {track}: wrote {len(selected)}/{limit} -> [cyan]{path}[/]")
        for record in selected[:3]:
            console.print(f"    {record['name']:44s} iptm={record['boltz_iptm_mean']:.3f} "
                          f"ph={record['ph_switch_score']:.1f} "
                          f"xr={record['cross_reactivity']}")


if __name__ == "__main__":
    app()
