"""Does histidine engineering help binding? Paired parent-vs-variant test across designs.

On one design all four variants beat the parent, but with high per-sample variance
(0.586 / 0.655 / 0.192 for one variant). This pairs every variant against its OWN parent,
folded under identical conditions, so design-to-design differences cancel and the question
becomes a within-family comparison.

Reports the binder-to-target ipTM only -- global ipTM is inflated by the TNF trimer's own
internal interfaces (~0.93) and is useless for ranking binders.
"""
import argparse
import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path


def binder_iptm(pred_dir: Path) -> tuple[float, float, float]:
    """(best, mean, min_ipde) over diffusion samples, binder-to-target pairs only."""
    best, means, ipdes = [], [], []
    for f in sorted(pred_dir.glob("confidence*.json")):
        c = json.loads(f.read_text())
        pc = c.get("pair_chains_iptm", {})
        if not pc:
            continue
        ids = sorted(pc.keys(), key=lambda x: int(x))
        b, tg = ids[0], ids[1:]
        vals = [float(pc[b][t]) for t in tg] + [float(pc[t][b]) for t in tg]
        best.append(max(vals))
        means.append(sum(vals) / len(vals))
        ipdes.append(float(c.get("complex_ipde", float("nan"))))
    if not best:
        return float("nan"), float("nan"), float("nan")
    return max(best), st.mean(means), min(ipdes)


ap = argparse.ArgumentParser()
ap.add_argument("manifest", type=Path)
ap.add_argument("pred_root", type=Path)
ap.add_argument("--out", type=Path)
args = ap.parse_args()

man = list(csv.DictReader(args.manifest.open()))
rows = []
for r in man:
    d = args.pred_root / r["spec"]
    if not d.is_dir():
        continue
    b, m, ipde = binder_iptm(d)
    if b != b:
        continue
    rows.append(r | {"best": b, "mean": m, "ipde": ipde})

fams: dict[str, dict] = defaultdict(dict)
for r in rows:
    if r["kind"] == "parent":
        fams[r["design"]]["parent"] = r
    else:
        fams[r["design"]].setdefault("variants", []).append(r)

print(f"{'design':26s} {'parent_best':>11s} {'parent_mean':>11s} "
      f"{'bestvar_best':>12s} {'bestvar_mean':>12s} {'delta_mean':>10s}  best variant")
deltas, improved, total = [], 0, 0
for design, fam in sorted(fams.items()):
    p = fam.get("parent")
    vs = fam.get("variants", [])
    if not p or not vs:
        continue
    bv = max(vs, key=lambda r: r["mean"])
    d = bv["mean"] - p["mean"]
    deltas.append(d)
    total += 1
    improved += d > 0
    print(f"{design:26s} {p['best']:11.3f} {p['mean']:11.3f} "
          f"{bv['best']:12.3f} {bv['mean']:12.3f} {d:+10.3f}  {bv['positions']}")

if deltas:
    print(f"\nPAIRED RESULT across {total} families:")
    print(f"  variants improved parent mean ipTM in {improved}/{total} families")
    print(f"  mean delta = {st.mean(deltas):+.3f}   median = {st.median(deltas):+.3f}")
    if len(deltas) > 1:
        print(f"  stdev of delta = {st.stdev(deltas):.3f}")
    pos = sum(1 for d in deltas if d > 0)
    print(f"  sign test: {pos} up / {len(deltas) - pos} down")

print("\nTop 12 constructs overall by binder-to-target ipTM (mean over samples):")
for r in sorted(rows, key=lambda x: -x["mean"])[:12]:
    tag = "parent " if r["kind"] == "parent" else "variant"
    print(f"  {r['spec']:46s} {tag} best={r['best']:.3f} mean={r['mean']:.3f} "
          f"ipde={r['ipde']:.2f}  {r['pkas']}")

if args.out:
    keys = ["spec", "design", "kind", "positions", "pkas", "best", "mean", "ipde", "sequence"]
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda x: -x["mean"]))
    print(f"\nwrote {args.out}")
