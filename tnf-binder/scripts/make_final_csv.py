"""Build the single ranked submission CSV in Proteinbase's actual round-2 format.

Corrections against the live submit page (the FAQ text differed):
  - molecule_class must be one of single_chain, nanobody, scfv, fab_kappa, fab_lambda.
    "protein" is NOT accepted; these de novo miniproteins are single_chain.
  - Up to 20 proteins per SUBMISSION, not 40/20/20 split across tracks.
  - 10-250 aa per chain.
  - Novelty is auto-scored after upload and must reach 3/4, so memorised sequences are
    actively harmful, not merely wasted slots.

Ranking is a single combined order, top row first, since there is one file rather than one per
track: pH switch first (the top-ranked objective, and the FAQ says a weak-but-selective binder
beats a strong non-selective one), then interface confidence, then cross-reactivity.
"""
import argparse
import csv
from pathlib import Path

VALID = {"single_chain", "nanobody", "scfv", "fab_kappa", "fab_lambda"}
MOLECULE_CLASS = "single_chain"
assert MOLECULE_CLASS in VALID
CAP = 20


def fnum(r, k, d=0.0):
    try:
        return float(r.get(k) or d)
    except (TypeError, ValueError):
        return d


ap = argparse.ArgumentParser()
ap.add_argument("tracks", type=Path, nargs="+", help="per-track CSVs to merge")
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()

seen, rows = set(), []
for t in args.tracks:
    for r in csv.DictReader(t.open()):
        s = r["sequence"].strip().upper()
        if s in seen:
            continue
        seen.add(s)
        rows.append(r)

rows.sort(key=lambda r: (-fnum(r, "ph_switch_score"), -fnum(r, "boltz_iptm_mean"),
                         -fnum(r, "cross_reactivity")))
rows = rows[:CAP]

# LF line endings to match the official template exactly. csv.writer defaults to CRLF
# (RFC 4180); the template is LF-only and their parser is untested here, so match it.
with args.out.open("w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["name", "sequence", "molecule_class"])
    for i, r in enumerate(rows, 1):
        seq = r["sequence"].strip().upper()
        assert 10 <= len(seq) <= 250, f"length {len(seq)} out of range"
        w.writerow([f"tnf_{i:02d}", seq, MOLECULE_CLASS])

print(f"wrote {len(rows)} designs -> {args.out} (cap {CAP})")
for i, r in enumerate(rows, 1):
    print(f"  {i:2d}. len={len(r['sequence']):3d} ph={r['ph_switch_score']:>4} "
          f"iptm={r['boltz_iptm_mean']:>6} xr={r['cross_reactivity']:>7} "
          f"his={r['pka_detail'] or '-'}")
