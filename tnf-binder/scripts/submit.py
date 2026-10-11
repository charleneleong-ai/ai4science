"""From every measurement to the one file that gets uploaded.

    submit.py shortlist --campaign scored.csv --outdir shortlists/
    submit.py rank --robustness pka.csv --human human/ --mouse mouse/ --out ranked.csv
    submit.py final shortlists/track2.csv --out submission.csv

Two stages because they answer different questions. The first applies the hard gates and ranks
what survives -- gates that DROP rather than down-rank, because the top-scoring construct of one
campaign was verbatim ubiquitin and the top of a 40-design shortlist was 59% alanine, and no
ranking rescues either. The second takes a shortlist and produces the capped, format-checked
file the competition accepts.

Both write against the rules in `tnf_binder.submission`, which is the same module
`tnf-binder validate` checks an upload with -- so a file this produces cannot fail a
pre-flight that this did not already catch.

Ranking follows the competition's own order: pH-selectivity, then mouse cross-reactivity, then
affinity, and its statement that a weak but clearly pH-selective binder beats a high-affinity
one that is not selective.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.confidence import DesignConfidence
from tnf_binder.design import TnfDesign
from tnf_binder.submission import (MAX_ROWS, MOLECULE_CLASS, SUBMISSION_COLUMNS,
                                   candidate_problem, row_problem, write_submission)
from tnf_binder.switch import SwitchProfile, switch_profiles
from tnf_binder.tables import design_stem, number, read_rows, write_rows
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

FILL_COLUMNS = ["design", "delta_q_floor", "his_robust", "human_iptm", "mouse_iptm",
                "mouse_retained", "interface_area", "cross_reactivity", "clashes",
                "flag", "sequence"]


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
    return f"structural flag: {flag}" if flag else None


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
def shortlist(
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


def missing_input(design: str, human: DesignConfidence, mouse: DesignConfidence,
                  geometry: dict[str, str] | None) -> str | None:
    """Which measurement this candidate is missing, if any.

    Reported separately from the gates and never folded into them. A missing geometry row
    leaves `clashes` at 0 and a missing fold directory leaves ipTM at 0.0 -- so one failure
    reads as "clash-free" and the other as "0% mouse retention", both of which are verdicts
    about the design rather than about our own inputs. `DesignConfidence.__bool__` is false for
    a directory that does not exist, which is the only way to tell a typo in `--human` from a
    binder that genuinely did not fold.
    """
    if not human:
        return f"no human fold under F__{design}"
    if not mouse:
        return f"no mouse fold under FM__{design}"
    if geometry is None:
        return "no geometry row (run score_designs first)"
    return None


def fill_record(design: str, profile: SwitchProfile, human: DesignConfidence,
                mouse: DesignConfidence, geometry: dict[str, str]) -> dict[str, object]:
    """One candidate, measured on all three objectives.

    `mouse_retained` is a ratio of two predictions of the same binder against the two orthologs,
    which is the only way objective 2 was ever measured here -- conserved-epitope coverage was
    tried as a proxy and was wrong on the design it had promoted (0.79 coverage, 47% retained).
    """
    retained = (mouse.binder_to_target / human.binder_to_target
                if human.binder_to_target else 0.0)
    return {
        "design": design,
        "delta_q_floor": round(profile.delta_q_floor, 2),
        "his_robust": profile.robust_tag,
        "human_iptm": round(human.binder_to_target, 3),
        "mouse_iptm": round(mouse.binder_to_target, 3),
        "mouse_retained": round(retained, 3),
        "interface_area": number(geometry, "interface_area") if geometry else "",
        "cross_reactivity": number(geometry, "cross_reactivity") if geometry else "",
        "clashes": int(number(geometry, "clashes")) if geometry else 0,
        "flag": (geometry.get("flag") or "").strip().strip('"'),
        "sequence": (geometry.get("sequence") or "").strip().upper(),
    }


@app.command()
def rank(
    robustness: Path = typer.Option(..., "--robustness",
                                    help="long-format pKa-per-sample CSV"),
    human: Path = typer.Option(..., "--human", help="Boltz predictions dir, human complexes"),
    mouse: Path = typer.Option(..., "--mouse", help="Boltz predictions dir, mouse complexes"),
    out: Path = typer.Option(..., "--out"),
    geometry: Path | None = typer.Option(None, "--geometry", help="table from score_designs"),
) -> None:
    """Rank slot-filling candidates on delta-Q, robustness and measured mouse cross-reactivity.

    Reports every candidate and why each one was dropped, not just the survivors: "3 of 57
    survived" is the result and the 3 alone is marketing. The gates are the module constants
    above, fixed before these designs were folded.
    """
    profiles = {design_stem(p.design): p for p in switch_profiles(read_rows(robustness))}
    # score_designs writes one row per diffusion sample; the gates are defined on model_0.
    geometry_by_design = {design_stem(r["design"]).removesuffix("_model_0"): r
                          for r in read_rows(geometry)
                          if r["design"].endswith("_model_0")} if geometry else {}
    console.print(f"{len(profiles)} candidates with pKa measured across samples")

    records, dropped, unmeasured = [], [], []
    for design, profile in sorted(profiles.items()):
        human_fold = DesignConfidence.from_dir(human / f"F__{design}")
        mouse_fold = DesignConfidence.from_dir(mouse / f"FM__{design}")
        geometry_row = geometry_by_design.get(design)
        if gap := missing_input(design, human_fold, mouse_fold, geometry_row):
            unmeasured.append((design, gap))
            continue
        record = fill_record(design, profile, human_fold, mouse_fold, geometry_row)
        if reason := candidate_problem(record):
            dropped.append((design, reason))
        else:
            records.append(record)

    records.sort(key=lambda r: (-r["delta_q_floor"], -r["mouse_retained"]))
    console.print(f"gates: kept [bold]{len(records)}[/] of {len(profiles)}, "
                  f"dropped [bold]{len(dropped)}[/], "
                  f"unmeasured [bold]{len(unmeasured)}[/]")
    for design, gap in unmeasured:
        console.print(f"  [yellow]{design}: {gap}[/]")
    for design, reason in dropped:
        console.print(f"  [dim]{design}: {reason}[/]")
    for record in records:
        console.print(f"  [green]{record['design']:26s}[/] dQ={record['delta_q_floor']:.2f} "
                      f"his={record['his_robust']} mouse={record['mouse_retained']:.0%} "
                      f"iptm={record['human_iptm']:.3f}")

    out.parent.mkdir(parents=True, exist_ok=True)
    write_rows(out, records, FILL_COLUMNS)
    console.print(f"wrote [bold]{len(records)}[/] ranked candidates -> [cyan]{out}[/]")


@app.command()
def final(
    tracks: list[Path] = typer.Argument(..., help="per-track CSVs to merge"),
    out: Path = typer.Option(..., "--out", help="the file that gets uploaded"),
    cap: int = typer.Option(MAX_ROWS, "--cap", help="rows per submission"),
) -> None:
    """Merge, rank and cap the shortlists into one upload-ready CSV."""
    seen, rows = set(), []
    for track in tracks:
        for row in read_rows(track):
            sequence = row["sequence"].strip().upper()
            if sequence in seen:
                continue
            seen.add(sequence)
            rows.append(row)

    rows.sort(key=lambda r: (-number(r, "ph_switch_score"), -number(r, "boltz_iptm_mean"),
                             -number(r, "cross_reactivity")))
    rows = rows[:cap]

    named = [(f"tnf_{i:02d}", r["sequence"].strip().upper()) for i, r in enumerate(rows, 1)]
    for name, sequence in named:
        if problem := row_problem({"name": name, "sequence": sequence,
                                   "molecule_class": MOLECULE_CLASS}):
            typer.echo(f"refusing to write: {problem}", err=True)
            raise typer.Exit(code=1)

    out.parent.mkdir(parents=True, exist_ok=True)
    write_submission(out, SUBMISSION_COLUMNS,
                     [[name, sequence, MOLECULE_CLASS] for name, sequence in named])
    console.print(f"wrote [bold]{len(named)}[/] designs -> [cyan]{out}[/] (cap {cap})")
    for (name, sequence), row in zip(named, rows):
        console.print(f"  {name} len={len(sequence):3d} ph={row['ph_switch_score']:>4} "
                      f"iptm={row['boltz_iptm_mean']:>6} xr={row['cross_reactivity']:>7} "
                      f"his={row['pka_detail'] or '-'}")


if __name__ == "__main__":
    app()
