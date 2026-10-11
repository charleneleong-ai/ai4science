"""Stage designs for the folding step.

WARNING -- do not rank by geometry here. Measured on 6 folded de novo designs, pre-fold
geometry does not predict interface confidence and may be anti-correlated with it:

  design  pre-fold geometry                      design_to_target_iptm
  msa_5   flagged: 171 clashes, low complexity   0.1797  (best)
  msa_1   clean, 414 A^2, epitope_focus 0.31     0.1615
  msa_3   cleanest, 604 A^2, focus 0.41          0.1145  (5th of 6)
  msa_2   ubiquitin, tiny interface              0.1086  (worst)

An earlier version of this file selected the top designs by epitope_focus, which would have
systematically folded the least promising ones. Geometry is still worth computing -- it feeds
cross-reactivity and the pH-switch scores, and it gates obvious junk out of SUBMISSIONS -- but
it is not a predictor of whether the model believes the interface.

So staging is deliberately unbiased: drop only designs that cannot be submitted at all
(memorised natural sequences, which are disqualified regardless of how well they fold), then
take an evenly spaced sample so a partial fold run still covers the batch.
"""
import argparse
import csv
import shutil
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("source_dir", type=Path, help="dir of .cif/.npz pairs")
    ap.add_argument("dest_dir", type=Path)
    ap.add_argument("--novelty", type=Path,
                    help="novelty_check CSV; memorised designs are dropped (disqualified)")
    ap.add_argument("--top", type=int, default=200, help="how many to stage")
    args = ap.parse_args()

    cifs = sorted(p for p in args.source_dir.rglob("*.cif") if "refold" not in p.parts)
    if not cifs:
        print(f"no designs under {args.source_dir}")
        return

    drop: set[str] = set()
    if args.novelty and args.novelty.exists():
        for r in csv.DictReader(args.novelty.open()):
            if r.get("novel", "").strip().lower() != "true":
                drop.add(r["design"])

    def memorised(stem: str) -> bool:
        return any(stem in d or d in stem for d in drop) if drop else False

    cand = [p for p in cifs if not memorised(p.stem)]
    if len(cand) > args.top:
        step = len(cand) / args.top
        cand = [cand[int(i * step)] for i in range(args.top)]

    args.dest_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for cif in cand:
        shutil.copy2(cif, args.dest_dir / cif.name)
        npz = cif.with_suffix(".npz")
        if npz.exists():
            shutil.copy2(npz, args.dest_dir / npz.name)
        n += 1
    print(f"{len(cifs)} designs | {len(drop)} memorised dropped | staged {n} (unbiased sample) "
          f"-> {args.dest_dir}")


if __name__ == "__main__":
    main()
