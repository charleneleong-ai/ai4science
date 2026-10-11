"""Mutate residues to histidine with a real, repacked side chain, via PyRosetta.

Needed because propka cannot score an incomplete residue: renaming a residue to HIS and keeping
only backbone+CB leaves no imidazole (ND1/CD2/CE1/NE2), so propka silently omits it and every
pKa comes back missing. PyRosetta builds the side chain and repacks the local neighbourhood, so
the predicted pKa reflects a physically plausible rotamer in its actual environment.

Runs in the BindCraft conda env, which is where pyrosetta lives. PYTHONNOUSERSITE=1 is required:
~/.local numpy shadows that env otherwise.
"""
import argparse
import sys
from pathlib import Path

import pyrosetta as pr
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb", type=Path)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--resnum", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--repack_radius", type=float, default=6.0)
    args = ap.parse_args()

    pr.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false", silent=True)
    pose = pr.pose_from_pdb(str(args.pdb))

    idx = pose.pdb_info().pdb2pose(args.chain, args.resnum)
    if idx == 0:
        print(f"ERROR: {args.chain}{args.resnum} not found in pose", file=sys.stderr)
        sys.exit(2)

    pr.toolbox.mutants.mutate_residue(pose, idx, "H", args.repack_radius,
                                      pr.get_score_function())
    pose.dump_pdb(str(args.out))
    name3 = pose.residue(idx).name3()
    print(f"OK {args.chain}{args.resnum} -> {name3} ; wrote {args.out}")
    if name3 != "HIS":
        print(f"WARNING: residue is {name3}, not HIS", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
