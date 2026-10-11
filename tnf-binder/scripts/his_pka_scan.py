"""Scan every interface position for histidine substitutions with a switchable pKa.

This replaces mechanism-based reasoning with a measurement. Objective 3 needs binding at pH 7.4
and release at pH 6.0, which requires a histidine whose pKa lies BETWEEN them: neutral at 7.4,
cationic at 6.0. Both mechanisms reasoned about earlier push pKa the WRONG way --

  PROPKA on a real design put three interface histidines at pKa 3.15, 4.70 and 5.71, all far
  below the 6.0 trigger. Interface burial destabilises the protonated form, and so does an
  adjacent cation. A buried His near a lysine is therefore the LEAST switchable position, yet
  it scored highest under the earlier structural proxy.

So: substitute His at every interface position in turn, predict pKa for each, and keep only
those landing in the usable window. PROPKA is CPU-only, so this runs without competing for GPU.
"""
import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

sys.path.insert(0, str(Path(__file__).parent))
from score_designs import load, classify_chains, residue_sasa, chain_seq
from pka_switch import run_propka, classify, IDEAL_LO, IDEAL_HI, USABLE_LO, USABLE_HI

MIN_INTERFACE_DSASA = 15.0
KEEP = set("GPC")          # Gly/Pro set backbone geometry, Cys may be bonded
_PYROSETTA_READY = False


def interface_positions(model, binders, targets) -> list[tuple]:
    """(chain_id, resnum, wt_aa, dsasa) for binder residues buried at the interface."""
    bids = {c.id for c in binders}
    tids = {c.id for c in targets}
    bound = residue_sasa(model, tids | bids)
    free = residue_sasa(model, bids)
    out = []
    for ch in binders:
        for r in (x for x in ch if x.id[0] == " "):
            key = (ch.id, r.id[1])
            d = free.get(key, 0.0) - bound.get(key, 0.0)
            aa = THREE2ONE.get(r.get_resname(), "X")
            if d >= MIN_INTERFACE_DSASA and aa not in KEEP and aa not in ("H", "X"):
                out.append((ch.id, r.id[1], aa, round(d, 1)))
    out.sort(key=lambda t: -t[3])
    return out


def write_his_mutant(src: Path, chain: str, resnum: int, dest: Path) -> bool:
    """Mutate one residue to HIS with a real, repacked side chain, via PyRosetta.

    Renaming the residue and keeping only backbone+CB is not enough: propka cannot score a
    histidine with no imidazole (ND1/CD2/CE1/NE2) and silently omits it, so every pKa came
    back missing. PyRosetta builds the side chain and repacks the neighbourhood.
    """
    import pyrosetta as pr

    global _PYROSETTA_READY
    if not _PYROSETTA_READY:
        pr.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false", silent=True)
        _PYROSETTA_READY = True
    pose = pr.pose_from_pdb(str(src))
    idx = pose.pdb_info().pdb2pose(chain, resnum)
    if idx == 0:
        return False
    pr.toolbox.mutants.mutate_residue(pose, idx, "H", 6.0, pr.get_score_function())
    if pose.residue(idx).name3() != "HIS":
        return False
    pose.dump_pdb(str(dest))
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", type=Path, nargs="+", help=".pdb complexes")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max_positions", type=int, default=14)
    args = ap.parse_args()

    paths: list[Path] = []
    for s in args.structures:
        paths.extend(sorted(p for p in s.rglob("*.pdb")) if s.is_dir() else [s])
    if not paths:
        print("no .pdb inputs (propka needs PDB)")
        return

    rows = []
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        for src in paths:
            model = load(src)
            targets, binders = classify_chains(model)
            if not binders or not targets:
                print(f"{src.stem}: chain classification failed")
                continue
            positions = interface_positions(model, binders, targets)[:args.max_positions]
            print(f"\n{src.stem}: {len(positions)} interface positions to scan")
            for ch, num, wt, dsasa in positions:
                mut = wd / f"{src.stem}_{ch}{num}H.pdb"
                if not write_his_mutant(src, ch, num, mut):
                    continue
                pkas = run_propka(mut, wd)
                pka = pkas.get((ch, num))
                if pka is None:
                    print(f"   {ch}{num} {wt}->H  propka gave no pKa")
                    continue
                verdict = classify(pka)
                rows.append({
                    "parent": src.stem, "chain": ch, "resnum": num, "wt": wt,
                    "dsasa": dsasa, "pka": round(pka, 2), "verdict": verdict,
                    "switchable": verdict in ("ideal", "usable"),
                })
                mark = {"ideal": "<== IDEAL", "usable": "<-- usable", "dead": ""}[verdict]
                print(f"   {ch}{num} {wt}->H  dSASA={dsasa:5.1f}  pKa={pka:5.2f}  "
                      f"{verdict:6s} {mark}")

    if not rows:
        print("\nno pKa results")
        return
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_ideal = sum(1 for r in rows if r["verdict"] == "ideal")
    n_use = sum(1 for r in rows if r["verdict"] == "usable")
    print(f"\nscanned {len(rows)} substitutions: {n_ideal} ideal "
          f"({IDEAL_LO}-{IDEAL_HI}), {n_use} usable ({USABLE_LO}-{USABLE_HI}) -> {args.out}")


if __name__ == "__main__":
    main()
