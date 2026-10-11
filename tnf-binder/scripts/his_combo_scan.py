"""Build multi-histidine pH switches and re-predict pKa for the COMBINATION.

Single-mutant pKa values do not compose. Two histidines placed near each other shift one
another's pKa (a protonated imidazole destabilises its neighbour's protonation), so a set of
positions that each look ideal alone can collapse once combined. Multiple switchable His are
still worth pursuing -- they give a steeper, cooperative transition between pH 7.4 and 6.0 than
any single residue -- but the combination has to be measured, not inferred.

For each combination this reports every introduced histidine's pKa in the mutant context, and
counts how many remain inside the switch window.
"""
import argparse
import csv
import itertools
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pka_switch import run_propka, classify
from score_designs import load, classify_chains, residue_sasa, chain_seq
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

_READY = False


def mutate_many(src: Path, muts: list[tuple[str, int]], dest: Path) -> bool:
    """Mutate several positions to HIS with repacked side chains."""
    import pyrosetta as pr

    global _READY
    if not _READY:
        pr.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false", silent=True)
        _READY = True
    pose = pr.pose_from_pdb(str(src))
    sfx = pr.get_score_function()
    for chain, num in muts:
        idx = pose.pdb_info().pdb2pose(chain, num)
        if idx == 0:
            return False
        pr.toolbox.mutants.mutate_residue(pose, idx, "H", 6.0, sfx)
        if pose.residue(idx).name3() != "HIS":
            return False
    pose.dump_pdb(str(dest))
    return True


def binder_sequence(pdb: Path) -> str:
    model = load(pdb)
    _, binders = classify_chains(model)
    return "".join(chain_seq(c) for c in binders)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("parent", type=Path)
    ap.add_argument("--positions", required=True,
                    help="comma list like A2,A5,A8,A15 (chain letter + resnum)")
    ap.add_argument("--max_his", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep_dir", type=Path, help="keep mutant PDBs here for refolding")
    args = ap.parse_args()

    pos = []
    for tok in args.positions.split(","):
        tok = tok.strip()
        pos.append((tok[0], int(tok[1:])))

    combos: list[tuple] = []
    for k in range(1, args.max_his + 1):
        combos += list(itertools.combinations(pos, k))
    print(f"{len(combos)} combinations from {len(pos)} positions")

    if args.keep_dir:
        args.keep_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        for combo in combos:
            tag = "_".join(f"{c}{n}H" for c, n in combo)
            mut = (args.keep_dir or wd) / f"{args.parent.stem}__{tag}.pdb"
            if not mutate_many(args.parent, list(combo), mut):
                print(f"  {tag}: mutation failed")
                continue
            pkas = run_propka(mut, wd)
            vals = {f"{c}{n}": pkas.get((c, n)) for c, n in combo}
            got = [v for v in vals.values() if v is not None]
            in_window = [v for v in got if classify(v) in ("ideal", "usable")]
            n_ideal = sum(1 for v in got if classify(v) == "ideal")
            rows.append({
                "design": f"{args.parent.stem}__{tag}",
                "parent": args.parent.stem,
                "n_his": len(combo),
                "positions": " ".join(f"{c}{n}" for c, n in combo),
                "pkas": " ".join(f"{k}={v:.2f}" if v is not None else f"{k}=NA"
                                 for k, v in vals.items()),
                "n_in_window": len(in_window),
                "n_ideal": n_ideal,
                "sequence": binder_sequence(mut),
            })
            print(f"  {tag:28s} pKa: {rows[-1]['pkas']:48s} "
                  f"in_window={len(in_window)}/{len(combo)} ideal={n_ideal}")

    if not rows:
        print("no results")
        return
    rows.sort(key=lambda r: (-r["n_ideal"], -r["n_in_window"], -r["n_his"]))
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    best = rows[0]
    print(f"\nbest: {best['design']} -- {best['n_ideal']} ideal / "
          f"{best['n_in_window']} in-window of {best['n_his']} His")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
