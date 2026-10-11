"""Assemble ranked submission CSVs from every measurement, applying the hard gates.

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

import argparse
import csv
from collections import Counter
from pathlib import Path

LIMITS = {1: 40, 2: 20, 3: 20}
MIN_LEN, MAX_LEN = 10, 250
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Proteinbase accepts only these; "protein" is NOT one of them.
MOLECULE_CLASS = "single_chain"
VALID_CLASSES = {"single_chain", "nanobody", "scfv", "fab_kappa", "fab_lambda"}

# Composition gate, independent of any model score.
MAX_SINGLE_AA_FRAC = 0.30
MIN_DISTINCT_AA = 12
MAX_AG_FRAC = 0.45

NOTE = ("epitope=conserved_TNFR2_site; pH switch = PROPKA-verified His pKa in 6.0-7.6; "
        "iptm is Boltz-2 binder-to-target (BoltzGen folding gave ~0.2 on the same designs "
        "-- models disagree, both reported)")


def fnum(row: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key) or default)
    except (TypeError, ValueError):
        return default


def load(path: Path | None) -> list[dict]:
    return list(csv.DictReader(path.open())) if path and path.exists() else []


def composition_problem(seq: str) -> str | None:
    """Degenerate sequences that no interface score should be trusted on."""
    if not seq:
        return "empty"
    counts = Counter(seq)
    aa, n = counts.most_common(1)[0]
    if n / len(seq) >= MAX_SINGLE_AA_FRAC:
        return f"low_complexity({aa}={n / len(seq):.0%})"
    if len(counts) < MIN_DISTINCT_AA:
        return f"only_{len(counts)}_distinct_aa"
    if (counts["A"] + counts["G"]) / len(seq) > MAX_AG_FRAC:
        return f"ala_gly_rich({(counts['A'] + counts['G']) / len(seq):.0%})"
    return None


def gate_reason(seq: str, novel: bool | None, geom_row: dict) -> str | None:
    """Why a construct is disqualified, or None if it survives every hard gate."""
    if not (MIN_LEN <= len(seq) <= MAX_LEN):
        return f"length {len(seq)}"
    if set(seq) - STD_AA:
        return "non-standard residues"
    if novel is False:
        return "MEMORISED (novelty gate)"
    if novel is None:
        return "no novelty verdict"
    problem = composition_problem(seq)
    if problem:
        return problem
    flag = (geom_row.get("flag") or "").strip()
    return f"structural flag [{flag}]" if flag else None


def ph_switch_terms(row: dict, pka_by_pos: dict) -> tuple[int, int, list[str]]:
    """PROPKA-verified in-window histidines for one construct."""
    n_ideal = n_usable = 0
    detail = []
    for token in (row.get("positions") or "").split():
        rec = pka_by_pos.get((row["design"], token))
        if not rec:
            continue
        detail.append(f"{token}:pKa={rec['pka']}")
        if rec["verdict"] == "ideal":
            n_ideal += 1
        elif rec["verdict"] == "usable":
            n_usable += 1
    return n_ideal, n_usable, detail


def novelty_lookup(novelty_rows: list[dict]):
    """spec -> bool, tolerant of the run/chain suffixes the various tools append."""
    if not novelty_rows:
        return lambda _spec: True
    verdict = {r["design"]: r["novel"].strip().lower() == "true" for r in novelty_rows}

    def check(spec: str) -> bool | None:
        if spec in verdict:
            return verdict[spec]
        return next((v for k, v in verdict.items() if spec in k or k in spec), None)

    return check


def build_record(row: dict, pka_by_pos: dict, geom_row: dict) -> dict:
    n_ideal, n_usable, detail = ph_switch_terms(row, pka_by_pos)
    return {
        "name": row["spec"],
        "sequence": (row.get("sequence") or "").strip().upper(),
        "boltz_iptm_best": round(fnum(row, "best"), 3),
        "boltz_iptm_mean": round(fnum(row, "mean"), 3),
        "boltz_ipde": round(fnum(row, "ipde"), 2),
        "n_his_ideal": n_ideal,
        "ph_switch_score": round(n_ideal + 0.4 * n_usable, 2),
        "pka_detail": ";".join(detail),
        "cross_reactivity": round(fnum(geom_row, "cross_reactivity"), 4) if geom_row else "",
        "epitope_focus": round(fnum(geom_row, "epitope_focus"), 4) if geom_row else "",
    }


def dedupe(records: list[dict]) -> list[dict]:
    seen, unique = set(), []
    for r in sorted(records, key=lambda x: -x["boltz_iptm_mean"]):
        if r["sequence"] not in seen:
            seen.add(r["sequence"])
            unique.append(r)
    return unique


def track_order(track: int):
    orders = {
        1: lambda r: (-r["boltz_iptm_mean"], -r["ph_switch_score"]),
        2: lambda r: (-(r["cross_reactivity"] or 0), -r["boltz_iptm_mean"]),
        3: lambda r: (-r["ph_switch_score"], -r["boltz_iptm_mean"]),
    }
    return orders[track]


def write_track(path: Path, records: list[dict], track: int, prefix: str) -> None:
    """LF line endings, matching the official template exactly."""
    columns = ["name", "sequence", "molecule_class", "boltz_iptm_best", "boltz_iptm_mean",
               "boltz_ipde", "n_his_ideal", "ph_switch_score", "pka_detail",
               "cross_reactivity", "epitope_focus", "notes"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        for i, r in enumerate(records, 1):
            writer.writerow([
                f"{prefix}_t{track}_{i:02d}", r["sequence"], MOLECULE_CLASS,
                r["boltz_iptm_best"], r["boltz_iptm_mean"], r["boltz_ipde"],
                r["n_his_ideal"], r["ph_switch_score"], r["pka_detail"],
                r["cross_reactivity"], r["epitope_focus"], f"src={r['name']}; {NOTE}",
            ])


def main() -> None:
    assert MOLECULE_CLASS in VALID_CLASSES
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--novelty", type=Path)
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--pka", type=Path, nargs="*", default=[])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="tnf")
    args = parser.parse_args()

    rows = load(args.campaign)
    pka_by_pos = {(r["parent"], f"{r['chain']}{r['resnum']}"): r
                  for path in args.pka for r in load(path)}
    geom = {r["design"]: r for r in load(args.geometry)}
    is_novel = novelty_lookup(load(args.novelty))
    print(f"{len(rows)} constructs from campaign")

    kept, dropped = [], []
    for row in rows:
        seq = (row.get("sequence") or "").strip().upper()
        geom_row = geom.get(row["spec"], {})
        reason = gate_reason(seq, is_novel(row["spec"]), geom_row)
        if reason:
            dropped.append(f"{row['spec']}: {reason}")
        else:
            kept.append(build_record(row, pka_by_pos, geom_row))

    print(f"gates: kept {len(kept)}, dropped {len(dropped)}")
    for message in dropped[:12]:
        print(f"  {message}")

    unique = dedupe(kept)
    print(f"unique sequences: {len(unique)}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    for track, limit in LIMITS.items():
        selected = sorted(unique, key=track_order(track))[:limit]
        path = args.outdir / f"{args.prefix}_track{track}.csv"
        write_track(path, selected, track, args.prefix)
        print(f"track {track}: wrote {len(selected)}/{limit} -> {path}")
        for r in selected[:3]:
            print(f"    {r['name']:44s} iptm={r['boltz_iptm_mean']:.3f} "
                  f"ph={r['ph_switch_score']:.1f} xr={r['cross_reactivity']}")


if __name__ == "__main__":
    main()
