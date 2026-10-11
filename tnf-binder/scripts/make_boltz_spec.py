"""Build a Boltz co-folding spec for a binder sequence against TNF-alpha.

This is the INDEPENDENT check on BindCraft. BindCraft optimises AF2's own interface ipTM, so
its reported 0.79 is self-confirming -- AF2 hallucination is known to be over-confident about
designs produced by backpropagating through AF2. Boltz is a different model and was not part of
the optimisation loop, so its verdict is the one that can be compared against the BoltzGen
numbers every other arm was measured with.

The target is given as the physiological TRIMER: three copies of the soluble domain. The
A+B dimer truncation measurably depressed interface confidence (target_ptm 0.47 vs 0.73).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq", required=True, help="binder sequence")
    ap.add_argument("--name", required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--protomers", type=int, default=3)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    chains = []
    chains.append(f"""  - protein:
      id: X
      sequence: {args.seq}
      msa: empty""")
    for cid in "ABC"[: args.protomers]:
        chains.append(f"""  - protein:
      id: {cid}
      sequence: {HUMAN}
      msa: empty""")
    spec = "version: 1\nsequences:\n" + "\n".join(chains) + "\n"
    p = args.out_dir / f"{args.name}.yaml"
    p.write_text(spec)
    print(f"wrote {p}  (binder {len(args.seq)} aa + {args.protomers} x TNF {len(HUMAN)} aa)")


if __name__ == "__main__":
    main()
