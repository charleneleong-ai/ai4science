"""Emit Boltz co-folding specs for a parent binder and its His variants.

Parent and variants are predicted the SAME way, so even though absolute binder-target ipTM is
low on this target (~0.2-0.3 across every method tried), the parent-vs-variant comparison is
still valid: it answers whether the engineered histidines destroy the interface.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN


def write_spec(out_dir: Path, name: str, seq: str, protomers: int = 3) -> None:
    chains = [f"""  - protein:
      id: X
      sequence: {seq}
      msa: empty"""]
    for cid in "ABC"[:protomers]:
        chains.append(f"""  - protein:
      id: {cid}
      sequence: {HUMAN}
      msa: empty""")
    (out_dir / f"{name}.yaml").write_text("version: 1\nsequences:\n" + "\n".join(chains) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("combo_csv", type=Path)
    ap.add_argument("--parent_seq", required=True)
    ap.add_argument("--parent_name", required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--pick", default="",
                    help="comma list of variant position-strings to include, e.g. 'A2 A5 A15'")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_spec(args.out_dir, args.parent_name, args.parent_seq)
    n = 1

    wanted = {w.strip() for w in args.pick.split(",") if w.strip()}
    for r in csv.DictReader(args.combo_csv.open()):
        if wanted and r["positions"] not in wanted:
            continue
        seq = r["sequence"].strip()
        if not seq or seq == args.parent_seq:
            continue
        tag = r["design"].split("__")[-1]
        write_spec(args.out_dir, f"var_{tag}", seq)
        n += 1
    print(f"wrote {n} specs (parent + variants) -> {args.out_dir}")
    for f in sorted(args.out_dir.glob("*.yaml")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
