"""Are the three diffusion samples of one design energetically distinguishable?

The folding model scores them at binder pTM 0.957 / 0.958 / 0.956 -- indistinguishable -- yet
PROPKA reads histidine 35 at pKa 3.63 on one and 9.69 on another. Either one placement is
strained and its pKa should be discounted, or the three are genuinely degenerate and three
draws badly undersample the conformational ensemble the switch depends on.

Rosetta REF15 answers the first half: it is parameterised for proteins, unlike the MLIP tier,
and is already installed. The second half is the local electrostatics, which is what PROPKA
actually responds to -- burial of the imidazole, and the charged side chains near it. A pKa of
3.6 means the neutral form is strongly favoured (desolvated, or beside a cation); 9.7 means the
protonated form is (beside a carboxylate).

Scores raw coordinates. No relax: the question is whether the structures *as predicted* differ,
and a relax would move the side chains under examination.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Selection
from Bio.PDB.SASA import ShrakeRupley

ACIDIC, BASIC = {"ASP", "GLU"}, {"LYS", "ARG"}
NEIGHBOUR_CUTOFF = 6.0
HIS_POSITIONS = (21, 35, 67)


def binder_chain(model):
    """The designed chain: 91 residues against the target's 157 each."""
    return min(model, key=lambda c: len(list(c)))


def environment(model, chain_id: str, resnum: int) -> dict[str, object]:
    """What sits around one histidine's imidazole, and how exposed it is."""
    his = model[chain_id][resnum]
    ring = [a for a in his if a.get_name() in ("ND1", "NE2", "CE1", "CD2", "CG")]
    if not ring:
        return {}
    centre = np.mean([a.coord for a in ring], axis=0)

    acidic = basic = 0
    nearest_acid = nearest_base = 99.0
    for residue in Selection.unfold_entities(model, "R"):
        if residue is his or residue.id[0] != " ":
            continue
        name = residue.get_resname()
        if name not in ACIDIC and name not in BASIC:
            continue
        d = min(float(np.linalg.norm(a.coord - centre)) for a in residue)
        if d <= NEIGHBOUR_CUTOFF:
            if name in ACIDIC:
                acidic += 1
            else:
                basic += 1
        if name in ACIDIC:
            nearest_acid = min(nearest_acid, d)
        else:
            nearest_base = min(nearest_base, d)
    return {"acidic_within_6A": acidic, "basic_within_6A": basic,
            "nearest_acidic": round(nearest_acid, 2), "nearest_basic": round(nearest_base, 2)}


def main() -> None:
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
    scorefxn = pyrosetta.get_fa_scorefxn()  # REF15

    parser, io, sr = MMCIFParser(QUIET=True), PDBIO(), ShrakeRupley()
    rows = []
    for cif in sorted(Path(sys.argv[1]).glob("*_model_*.cif")):
        structure = parser.get_structure(cif.stem, str(cif))
        model = structure[0]
        chain = binder_chain(model).id

        pdb = Path("/tmp") / f"{cif.stem}.pdb"
        io.set_structure(structure)
        io.save(str(pdb))
        pose = pyrosetta.pose_from_file(str(pdb))
        total = scorefxn(pose)

        sr.compute(model, level="R")
        row = {"sample": cif.stem[-7:], "binder_chain": chain, "residues": pose.total_residue(),
               "ref15_total": round(total, 1),
               "ref15_per_res": round(total / pose.total_residue(), 3)}
        for p in HIS_POSITIONS:
            his = model[chain][p]
            row[f"his{p}_sasa"] = round(his.sasa, 1)
            row[f"his{p}_env"] = environment(model, chain, p)
        rows.append(row)

    print(f"{'sample':<9}{'chain':>6}{'res':>6}{'REF15 total':>13}{'per-res':>9}")
    for r in rows:
        print(f"{r['sample']:<9}{r['binder_chain']:>6}{r['residues']:>6}"
              f"{r['ref15_total']:>13.1f}{r['ref15_per_res']:>9.3f}")
    spread = max(r["ref15_total"] for r in rows) - min(r["ref15_total"] for r in rows)
    print(f"\nREF15 spread across the three samples: {spread:.1f} REU "
          f"({spread / rows[0]['residues']:.3f} per residue)")

    for p in HIS_POSITIONS:
        print(f"\nHis {p} — side-chain SASA and what is within {NEIGHBOUR_CUTOFF} A of the ring:")
        for r in rows:
            e = r[f"his{p}_env"]
            print(f"  {r['sample']}  SASA {r[f'his{p}_sasa']:>6.1f} A^2   "
                  f"acidic {e['acidic_within_6A']}  basic {e['basic_within_6A']}   "
                  f"nearest acid {e['nearest_acidic']:>5.2f} A  nearest base {e['nearest_basic']:>5.2f} A")


if __name__ == "__main__":
    main()
