"""Compare screen configs on design_to_target_iptm, the binder-specific interface score.

Reports per-config statistics so a winner is chosen on the whole batch rather than on a single
design -- an n=1 read earlier reversed itself at n=3.

Metrics are stored per diffusion sample (shape (5,)); this uses the BEST sample per design,
since that is the structure one would keep, and also shows the per-sample mean.

Deliberately does NOT rank on global iptm: that averages over all inter-chain pairs, so a
multi-chain target inflates it with its own well-predicted internal interfaces (dimer -> trimer
moved global iptm 0.20 -> 0.52 while the binder interface stayed flat).
"""
import statistics as st
import sys
from pathlib import Path

import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("out")
CONFIGS = sys.argv[2].split(",") if len(sys.argv) > 2 else [
    "screen_pocket_free", "screen_pocket_helix", "screen_pocket5_helix", "screen_lobe3_free",
]
# Reference arms already measured, for context.
BASELINES = {
    "nano (dimer)": ROOT / "nano/intermediate_designs_inverse_folded/fold_out_npz",
    "nano (trimer)": ROOT / "nano_trimer/intermediate_designs_inverse_folded/fold_out_npz",
    "denovo (dimer)": ROOT / "v2_msa/intermediate_designs_inverse_folded/fold_out_npz",
}


def collect(d: Path) -> list[dict[str, float]]:
    rows = []
    if not d or not d.exists():
        return rows
    for p in sorted(d.rglob("*.npz")):
        try:
            z = np.load(p, allow_pickle=True)
        except Exception:
            continue
        def g(k, how="max"):
            # Metrics are per diffusion sample (shape (5,)). "best" is the sample one would
            # keep; "mean" describes the config. Never take element 0 -- samples trend.
            if k not in z.files:
                return float("nan")
            a = np.asarray(z[k]).reshape(-1).astype(float)
            if not a.size:
                return float("nan")
            return float(a.max() if how == "max" else a.min() if how == "min" else a.mean())
        rows.append({"name": p.stem,
                     "d2t": g("design_to_target_iptm", "max"),
                     "d2t_mean": g("design_to_target_iptm", "mean"),
                     "iptm": g("iptm", "max"),
                     "pae": g("interaction_pae", "min"),
                     "target_ptm": g("target_ptm", "mean")})
    return rows


def summarise(label: str, rows: list[dict[str, float]]) -> None:
    vals = [r["d2t"] for r in rows if r["d2t"] == r["d2t"]]
    if not vals:
        print(f"  {label:24s} n=0   (no folded designs yet)")
        return
    best = max(rows, key=lambda r: (r["d2t"] if r["d2t"] == r["d2t"] else -1))
    paes = [r["pae"] for r in rows if r["pae"] == r["pae"]]
    means = [r["d2t_mean"] for r in rows if r["d2t_mean"] == r["d2t_mean"]]
    print(f"  {label:24s} n={len(vals):<3} d2t_best_mean={st.mean(vals):.3f}  "
          f"d2t_top={max(vals):.3f}  (per-config mean={st.mean(means):.3f})  "
          f"pae_min={min(paes):.1f}  top={best['name']}")


print("=== reference arms ===")
for label, d in BASELINES.items():
    summarise(label, collect(d))

print("\n=== screen configs (trimer target) ===")
results = {}
for c in CONFIGS:
    d = next((p for p in (ROOT / c).rglob("fold_out_npz")), None) if (ROOT / c).exists() else None
    rows = collect(d)
    results[c] = rows
    summarise(c, rows)

ranked = [(c, st.mean([r["d2t"] for r in rows if r["d2t"] == r["d2t"]]))
          for c, rows in results.items() if any(r["d2t"] == r["d2t"] for r in rows)]
if ranked:
    ranked.sort(key=lambda x: -x[1])
    print(f"\nwinner by mean design_to_target_iptm: {ranked[0][0]} ({ranked[0][1]:.3f})")
    if len(ranked) > 1:
        print(f"  runner-up: {ranked[1][0]} ({ranked[1][1]:.3f})")
    print("  bar for a credible de novo interface is ~0.5")
