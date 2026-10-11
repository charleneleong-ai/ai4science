"""Join BoltzGen folding confidence onto the geometric scores, and gate on both.

Two complementary views of a design, and both are needed:

  geometry (score_designs)  measures the pose BoltzGen DESIGNED -- interface area, how much of
                            it lands on the intended epitope, how much on species-divergent
                            residues.
  confidence (folding step) measures whether the model BELIEVES that pose when it re-predicts
                            the complex from scratch -- iptm, design_to_target_iptm,
                            interaction_pae.

A design can have excellent geometry and still be a non-binder: the first nanobody batch buried
plausible interfaces yet returned iptm ~0.15 with interaction_pae ~28 A. Confidence is the
stricter and more predictive filter, so it gates; geometry then orders within the survivors.

Gate on design_to_target_iptm, NOT the global iptm. Global iptm averages over every inter-chain
pair, so a multi-chain target inflates it with its own well-predicted internal interfaces: moving
from the TNF dimer to the trimer lifted global iptm 0.20 -> 0.52 purely by adding three
well-predicted TNF-TNF pairs, while design_to_target_iptm stayed flat at 0.19 and interaction_pae
did not improve. Gating on global iptm would have passed a design whose binder does not dock.

Thresholds follow de novo binder-design practice (BindCraft-style): >= 0.5 for the binder-to-target
interface and interaction pae <= 10. They are explicit and overridable rather than tuned to
whatever this run happened to produce.
"""
import argparse
import csv
from pathlib import Path


def read_csv(p: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(p.open())) if p and p.exists() else []


def fnum(r: dict[str, str], k: str, default: float = 0.0) -> float:
    try:
        return float(r.get(k) or default)
    except (TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scored_csv", type=Path)
    ap.add_argument("fold_csv", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min_d2t_iptm", type=float, default=0.5,
                    help="threshold on design_to_target_iptm, the binder-specific interface score")
    ap.add_argument("--max_pae", type=float, default=10.0)
    ap.add_argument("--relax", action="store_true",
                    help="keep designs failing the confidence gate, flagged, for inspection")
    args = ap.parse_args()

    scored = read_csv(args.scored_csv)
    fold = {r["design"]: r for r in read_csv(args.fold_csv)}

    out = []
    for r in scored:
        f = fold.get(r["design"], {})
        iptm = fnum(f, "iptm")
        d2t = fnum(f, "design_to_target_iptm")
        pae = fnum(f, "interaction_pae", 99.0)
        flags = [x for x in (r.get("flag") or "").split(",") if x]
        if not f:
            flags.append("no_fold_metrics")
        else:
            if d2t < args.min_d2t_iptm:
                flags.append(f"low_d2t_iptm({d2t:.2f})")
            if pae > args.max_pae:
                flags.append(f"high_pae({pae:.1f})")
        out.append(r | {
            "iptm": iptm, "design_to_target_iptm": d2t, "interaction_pae": pae,
            "target_ptm": fnum(f, "target_ptm"), "design_ptm": fnum(f, "design_ptm"),
            "flag": ",".join(flags),
        })

    if not args.relax:
        kept = [r for r in out if not r["flag"]]
    else:
        kept = out

    cols = ["design", "binder_len", "iptm", "design_to_target_iptm", "interaction_pae",
            "design_ptm", "target_ptm", "interface_area", "epitope_focus", "cross_reactivity",
            "divergent_frac", "n_his_interface", "ph_switch", "clashes", "flag", "sequence"]
    kept.sort(key=lambda r: (-fnum(r, "design_to_target_iptm"),
                             fnum(r, "interaction_pae", 99.0)))
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)

    n_pass = sum(1 for r in out if not r["flag"])
    print(f"{len(out)} scored | {n_pass} pass both gates "
          f"(design_to_target_iptm>={args.min_d2t_iptm}, pae<={args.max_pae}) | "
          f"wrote {len(kept)} -> {args.out}")
    if n_pass == 0 and out:
        best = max(out, key=lambda r: fnum(r, "design_to_target_iptm"))
        print(f"  best by design_to_target_iptm: {best['design']} "
              f"d2t={best['design_to_target_iptm']:.3f} (global iptm={best['iptm']:.3f}) "
              f"pae={best['interaction_pae']:.1f} -- nothing clears the bar")


if __name__ == "__main__":
    main()
