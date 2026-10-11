"""Emit Boltz specs for BindCraft designs so they can be scored on the same yardstick.

BindCraft reports AF2 i_pTM on designs it optimised against AF2, which is self-confirming: the
one design checked independently scored i_pTM 0.92 under AF2 but 0.21-0.26 under Boltz. The
campaign designs were all measured under Boltz, so these must be too before any comparison.
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


ap = argparse.ArgumentParser()
ap.add_argument("stats_csv", type=Path)
ap.add_argument("--out_dir", type=Path, required=True)
ap.add_argument("--min_iptm", type=float, default=0.6)
ap.add_argument("--manifest", type=Path, required=True)
args = ap.parse_args()

args.out_dir.mkdir(parents=True, exist_ok=True)
rows, kept = list(csv.DictReader(args.stats_csv.open())), []
for r in rows:
    seq = (r.get("Sequence") or "").strip().upper()
    try:
        af_iptm = float(r.get("i_pTM") or 0)
    except ValueError:
        af_iptm = 0.0
    if not seq or af_iptm < args.min_iptm:
        continue
    name = f"BC__{r['Design']}"
    write_spec(args.out_dir, name, seq)
    kept.append({"spec": name, "design": r["Design"], "kind": "bindcraft",
                 "af2_iptm": af_iptm, "af2_ipae": r.get("i_pAE"),
                 "shape_comp": r.get("ShapeComplementarity"),
                 "n_interface": r.get("n_InterfaceResidues"),
                 "rosetta_energy": r.get("Binder_Energy_Score"),
                 "positions": "", "pkas": "", "sequence": seq})

with args.manifest.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(kept[0].keys()) if kept else ["spec"])
    w.writeheader()
    w.writerows(kept)
print(f"{len(rows)} trajectories, {len(kept)} with AF2 i_pTM >= {args.min_iptm}")
for k in kept:
    print(f"  {k['spec']:30s} af2_iptm={k['af2_iptm']} len={len(k['sequence'])}")
print(f"wrote specs -> {args.out_dir}, manifest {args.manifest}")
