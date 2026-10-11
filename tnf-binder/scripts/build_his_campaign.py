"""Build the full His-switch campaign: every design with switchable positions, parent + variants.

Rationale from the measured data:
  - 200 single substitutions scanned, 18 ideal (pKa 6.3-7.2) and 17 usable (6.0-7.6).
  - No structural descriptor predicts which positions work (burial r=-0.37, carboxylate count
    r=+0.21, cation count the wrong sign), so positions come from the pKa scan, not a rule.
  - On the one design tested so far, ALL FOUR His variants improved binder-target ipTM over the
    parent (0.288 -> 0.366-0.444 mean). That is n=1, so this replicates it across designs.

Emits Boltz co-folding specs for each parent and its variants, so parent and variant are
always measured under identical conditions.
"""
import argparse
import csv
import itertools
import sys
from collections import defaultdict
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


def switchable_by_design(scan_paths) -> dict:
    """design -> switchable positions, ideal first, then nearest the middle of the window.

    Positions come from the pKa scan, never from a structural rule: burial correlates only
    r=-0.37 with predicted pKa, and carboxylate and cation counts fare worse still.
    """
    from collections import defaultdict
    found = defaultdict(list)
    for path in scan_paths:
        for r in csv.DictReader(path.open()):
            if r.get("switchable", "").strip().lower() == "true":
                found[r["parent"]].append(
                    (r["chain"], int(r["resnum"]), r["wt"], float(r["pka"]), r["verdict"]))
    for design in found:
        found[design].sort(key=lambda t: (t[4] != "ideal", abs(t[3] - 6.7)))
    return found


def binder_chain_of(pose, info) -> str:
    """Shortest chain is the binder; TNF protomers are 152 aa each."""
    from collections import defaultdict
    lengths = defaultdict(int)
    for i in range(1, pose.total_residue() + 1):
        lengths[info.chain(i)] += 1
    return min(lengths, key=lambda c: lengths[c])


def chain_sequence(pose, info, chain: str) -> str:
    return "".join(pose.residue(i).name1()
                   for i in range(1, pose.total_residue() + 1)
                   if info.chain(i) == chain)


def variant_combos(positions: list, max_his: int, cap: int) -> list:
    combos = []
    for k in range(1, min(max_his, len(positions)) + 1):
        combos += list(itertools.combinations(positions, k))
    return combos[:cap]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scans", type=Path, nargs="+", help="pka_scan csv files")
    ap.add_argument("--pdb_dirs", type=Path, nargs="+", required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--max_his", type=int, default=3)
    ap.add_argument("--max_variants_per_design", type=int, default=4)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()

    by_design = switchable_by_design(args.scans)

    print(f"{len(by_design)} designs carry switchable positions")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def find_pdb(stem: str) -> Path | None:
        for d in args.pdb_dirs:
            p = d / f"{stem}.pdb"
            if p.exists():
                return p
        return None

    import pyrosetta as pr
    pr.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false", silent=True)
    sfx = pr.get_score_function()

    rows = []
    for design, positions in sorted(by_design.items()):
        src = find_pdb(design)
        if src is None:
            print(f"  {design}: no pdb found, skipping")
            continue
        pose0 = pr.pose_from_pdb(str(src))
        info = pose0.pdb_info()
        binder_chain = binder_chain_of(pose0, info)
        parent_seq = chain_sequence(pose0, info, binder_chain)
        write_spec(args.out_dir, f"P__{design}", parent_seq)
        rows.append({"spec": f"P__{design}", "design": design, "kind": "parent",
                     "positions": "", "pkas": "", "sequence": parent_seq})

        usable = [p for p in positions if p[0] == binder_chain][:4]
        combos = variant_combos(usable, args.max_his, args.max_variants_per_design)

        for combo in combos:
            pose = pose0.clone()
            ok = True
            for ch, num, _wt, _pka, _v in combo:
                idx = info.pdb2pose(ch, num)
                if idx == 0:
                    ok = False
                    break
                pr.toolbox.mutants.mutate_residue(pose, idx, "H", 6.0, sfx)
            if not ok:
                continue
            tag = "_".join(f"{c}{n}H" for c, n, _, _, _ in combo)
            seq = chain_sequence(pose, info, binder_chain)
            name = f"V__{design}__{tag}"
            write_spec(args.out_dir, name, seq)
            rows.append({"spec": name, "design": design, "kind": "variant",
                         "positions": " ".join(f"{c}{n}" for c, n, _, _, _ in combo),
                         "pkas": " ".join(f"{c}{n}={pk:.2f}" for c, n, _, pk, _ in combo),
                         "sequence": seq})
        print(f"  {design}: binder chain {binder_chain}, "
              f"{len(usable)} switchable -> {len(combos)} variants")

    with args.manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["spec", "design", "kind", "positions", "pkas", "sequence"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} specs -> {args.out_dir}  (manifest {args.manifest})")


if __name__ == "__main__":
    main()
