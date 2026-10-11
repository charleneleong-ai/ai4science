"""Test what actually predicts a switchable histidine pKa.

Three candidate explanations, from the pKa scan data:
  burial          interface dSASA -- the thing the earlier structural proxy used
  carboxylate     nearby Asp/Glu carboxylate oxygens, which RAISE His pKa
  cation          nearby Lys/Arg nitrogens, which LOWER it

The scan showed deeply buried positions collapsing to pKa ~3 while moderately buried ones
reached 6.3-7.1, so burial alone does not explain it. This correlates each candidate against
the measured pKa to find the rule worth designing to.
"""
import csv
import sys
from pathlib import Path

import numpy as np
from Bio import PDB

sys.path.insert(0, str(Path(__file__).parent))

CARBOXYL = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATION = {("LYS", "NZ"), ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE")}
RADIUS = 8.0

rows = list(csv.DictReader(Path(sys.argv[1]).open()))
pdb_dir = Path(sys.argv[2])

out = []
cache: dict[str, object] = {}
for r in rows:
    parent = r["parent"]
    if parent not in cache:
        p = pdb_dir / f"{parent}.pdb"
        if not p.exists():
            continue
        cache[parent] = PDB.PDBParser(QUIET=True).get_structure(parent, p)[0]
    model = cache[parent]
    ch, num = r["chain"], int(r["resnum"])
    try:
        res = model[ch][num]
    except KeyError:
        continue
    cen = np.mean([a.coord for a in res if a.element != "H"], axis=0)

    n_cox = n_cat = 0
    for c in model:
        for rr in c:
            if rr.id[0] != " ":
                continue
            for a in rr:
                key = (rr.get_resname(), a.get_id())
                if key in CARBOXYL or key in CATION:
                    if np.linalg.norm(a.coord - cen) <= RADIUS:
                        if key in CARBOXYL:
                            n_cox += 1
                        else:
                            n_cat += 1
    out.append({
        "design": f"{parent}:{ch}{num}", "pka": float(r["pka"]),
        "dsasa": float(r["dsasa"]), "n_carboxyl": n_cox, "n_cation": n_cat,
        "net": n_cox - n_cat, "verdict": r["verdict"],
    })

if not out:
    print("no rows matched")
    sys.exit()

print(f"{'position':34s} {'pKa':>6s} {'dSASA':>7s} {'COO-':>5s} {'cat+':>5s} {'net':>5s} verdict")
for o in sorted(out, key=lambda x: -x["pka"]):
    print(f"{o['design']:34s} {o['pka']:6.2f} {o['dsasa']:7.1f} "
          f"{o['n_carboxyl']:5d} {o['n_cation']:5d} {o['net']:5d} {o['verdict']}")

def corr(key):
    x = np.array([o[key] for o in out], dtype=float)
    y = np.array([o["pka"] for o in out], dtype=float)
    if x.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])

print(f"\nPearson r vs predicted pKa  (n={len(out)}):")
for k in ("dsasa", "n_carboxyl", "n_cation", "net"):
    print(f"  {k:12s} r = {corr(k):+.3f}")
