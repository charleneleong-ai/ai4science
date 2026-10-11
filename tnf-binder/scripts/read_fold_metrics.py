"""Collect BoltzGen folding-step confidence metrics.

Each metric is stored per DIFFUSION SAMPLE, not as a scalar: coords has shape (5, n_atoms, 3)
and every confidence key is shape (5,). The samples are not interchangeable -- they trend
across the schedule, e.g. design_to_target_iptm [0.237, 0.287, 0.325, 0.352, 0.309] for one
design -- so collapsing them to a single number silently picks a statistic.

Both are reported:
  *_mean  describes the config, for comparing settings
  *_best  the sample you would actually keep, so it is what ranking uses

An earlier version returned the mean for size>1 and the first element elsewhere, which made two
readers disagree on the same file. Reporting both removes the ambiguity.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

# interface confidence: higher is better / error: lower is better
HIGHER = ("iptm", "design_iptm", "design_to_target_iptm", "design_iiptm",
          "ptm", "design_ptm", "target_ptm", "protein_iptm")
LOWER = ("interaction_pae", "min_interaction_pae", "min_design_to_target_pae")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("npz_dir", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--sort", default="design_to_target_iptm_best")
    args = ap.parse_args()

    paths = sorted(args.npz_dir.rglob("*.npz"))
    if not paths:
        print(f"no npz under {args.npz_dir}")
        return

    rows, nsamples = [], set()
    for p in paths:
        try:
            z = np.load(p, allow_pickle=True)
        except Exception as e:
            print(f"skip {p.name}: {e}", file=sys.stderr)
            continue
        row = {"design": p.stem}
        if "coords" in z.files:
            nsamples.add(int(np.asarray(z["coords"]).shape[0]))
        for k in HIGHER + LOWER:
            if k not in z.files:
                continue
            a = np.asarray(z[k]).reshape(-1).astype(float)
            if not a.size:
                continue
            row[f"{k}_mean"] = round(float(a.mean()), 4)
            row[f"{k}_best"] = round(float(a.max() if k in HIGHER else a.min()), 4)
        rows.append(row)

    if nsamples:
        print(f"diffusion samples per design: {sorted(nsamples)}", file=sys.stderr)

    cols = ["design"] + [c for c in
                         [f"{k}_{s}" for k in HIGHER + LOWER for s in ("best", "mean")]
                         if any(c in r for r in rows)]
    rows.sort(key=lambda r: -(r.get(args.sort) or -1e9))
    w = csv.DictWriter(sys.stdout, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
    if args.out:
        with args.out.open("w", newline="") as f:
            ww = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            ww.writeheader()
            ww.writerows(rows)
        print(f"\nwrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
