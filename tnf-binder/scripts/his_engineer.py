"""Propose histidine substitutions that make a binder release TNF-alpha at pH 6.0.

Objective 3 wants binding at pH 7.4 and none at pH 6.0: the NEUTRAL histidine must bind and
the CATIONIC one must not. Two mechanisms are worth engineering, and the intuitive
"His beside a lysine" idea is not one of them:

  buried  a His that is buried BY THE INTERFACE pays a desolvation penalty when it protonates
          in an apolar environment. This is the dominant mechanism in recycling antibodies.
  anchor  a His held 5-9 A from a conserved cationic anchor (R32/K90/K128) is repelled once
          protonated. But putting His in direct contact with a cation *lowers* its pKa, often
          to ~5.0-5.5, so the switch stops firing above pH 6. Contacts closer than
          MIN_ANCHOR_DIST are therefore rejected, not rewarded.

Burial is measured as interface dSASA (unbound minus bound), never as bound SASA. Bound SASA
counts the binder's own hydrophobic core as "buried", which produced 77 bogus candidate sites
on a 93-residue design in testing.

Clusters of 2-3 His raise the apparent pKa cooperatively, so multi-mutants are emitted too.
No folding model represents protonation, so these are structural proposals to be re-folded and
re-scored, never a free-energy claim.
"""
import argparse
import csv
import itertools
from pathlib import Path
import sys

import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

sys.path.insert(0, str(Path(__file__).parent))
from score_designs import load, classify_chains, chain_seq, mature_map, residue_sasa
from tnf_common import PH_ANCHORS

MIN_ANCHOR_DIST = 5.0      # closer and the His pKa collapses below the pH 6.0 trigger
MAX_ANCHOR_DIST = 9.0      # further and Coulombic coupling is negligible
MIN_INTERFACE_DSASA = 15.0  # A^2 a residue must bury at the interface to count as interfacial
KEEP = set("GPC")           # Gly/Pro shape the backbone, Cys may be in a disulfide
MAX_HIS = 3
ANCHOR_ATOMS = ("NZ", "NH1", "NH2", "NE")
HIS_REF_SASA = 180.0


def candidates(path: Path) -> tuple[str, list[dict]]:
    model = load(path)
    targets, binders = classify_chains(model)
    if not targets or not binders:
        return "", []

    tids = {c.id for c in targets}
    bids = {c.id for c in binders}
    bound = residue_sasa(model, tids | bids)
    b_free = residue_sasa(model, bids)

    tnum: dict = {}
    for ch in targets:
        tnum.update(mature_map(ch))
    anchors = [a for ch in targets for r in ch if r.id[0] == " "
               and tnum.get(r) in PH_ANCHORS for a in r if a.get_id() in ANCHOR_ATOMS]

    seq = "".join(chain_seq(ch) for ch in binders)
    out, idx = [], 0
    for ch in binders:
        for r in (x for x in ch if x.id[0] == " "):
            idx += 1
            aa = THREE2ONE.get(r.get_resname(), "X")
            if aa in KEEP or aa in ("H", "X"):
                continue
            key = (ch.id, r.id[1])
            dsasa = b_free.get(key, 0.0) - bound.get(key, 0.0)
            if dsasa < MIN_INTERFACE_DSASA:
                continue                        # not actually at the interface
            heavy = [a for a in r if a.element != "H"]
            d = min((np.linalg.norm(a.coord - b.coord) for a in heavy for b in anchors),
                    default=np.inf)
            mech = []
            if dsasa >= 0.4 * HIS_REF_SASA:
                mech.append("buried")
            if MIN_ANCHOR_DIST <= d <= MAX_ANCHOR_DIST:
                mech.append("anchor")
            if not mech:
                continue
            out.append({"pos": idx, "wt": aa, "dsasa": round(dsasa, 1),
                        "d_anchor": round(float(d), 2) if np.isfinite(d) else None,
                        "mechanism": "+".join(mech)})
    out.sort(key=lambda c: (-len(c["mechanism"].split("+")), -c["dsasa"]))
    return seq, out


def mutate(seq: str, positions: tuple[int, ...]) -> str:
    s = list(seq)
    for p in positions:
        s[p - 1] = "H"
    return "".join(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", type=Path)
    ap.add_argument("out_csv", type=Path)
    ap.add_argument("--max_per_design", type=int, default=8)
    ap.add_argument("--top_sites", type=int, default=4)
    args = ap.parse_args()

    rows = []
    for path in sorted(p for p in args.structures.rglob("*") if p.suffix in (".cif", ".pdb")):
        seq, cands = candidates(path)
        if not seq or not cands:
            print(f"{path.stem}: no viable His site")
            continue
        sites = cands[:args.top_sites]
        print(f"{path.stem}: {len(cands)} interface sites, using top {len(sites)}")
        for c in sites:
            print(f"   pos {c['pos']:>3} {c['wt']}->H  dSASA={c['dsasa']:>5.1f} "
                  f"d_anchor={c['d_anchor']} ({c['mechanism']})")
        by_pos = {c["pos"]: c for c in sites}
        combos: list[tuple[int, ...]] = []
        for k in range(1, min(MAX_HIS, len(sites)) + 1):
            combos += list(itertools.combinations([c["pos"] for c in sites], k))
        for combo in combos[:args.max_per_design]:
            rows.append({
                "design": f"{path.stem}_H{'_'.join(map(str, combo))}",
                "parent": path.stem,
                "sequence": mutate(seq, combo),
                "n_his_added": len(combo),
                "his_positions": " ".join(map(str, combo)),
                "mechanisms": ";".join(by_pos[p]["mechanism"] for p in combo),
                "mean_dsasa": round(float(np.mean([by_pos[p]["dsasa"] for p in combo])), 1),
            })

    if not rows:
        print("no mutants proposed")
        return
    with args.out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nproposed {len(rows)} His variants -> {args.out_csv}")


if __name__ == "__main__":
    main()
