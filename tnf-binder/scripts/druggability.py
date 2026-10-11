"""Rank candidate hotspot clusters on how designable they are, not just how conserved.

Conservation told us WHERE we are allowed to bind. It says nothing about whether a binder can
grip there. A shallow, polar, solvent-exposed groove is a hard de novo target; a concave
hydrophobic pocket is an easy one. Since every architecture tried so far plateaued at
design_to_target_iptm 0.11-0.30, the epitope's own designability is a prime suspect.

For each conserved, exposed epitope residue this scores a tight neighbourhood on:
  concavity    mean count of target heavy atoms within 10 A of the residue's outward normal
               probe point -- high means the probe sits in a pocket, low means a flat surface
  hydrophobic  fraction of neighbourhood SASA contributed by apolar side chains
  compactness  max pairwise span of the cluster (tight clusters are easier to satisfy)
"""
import sys
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import divergent_positions

SD = Path(sys.argv[1])
DIVERGENT = divergent_positions()
APOLAR = set("AVLIMFWCYP")
# All TNFR2 contact positions (mature numbering), from the 3ALQ analysis.
EPITOPE_ALL = [20, 21, 23, 31, 32, 33, 63, 67, 71, 72, 73, 75, 77, 85, 86, 87, 88, 89, 90,
               91, 92, 97, 111, 113, 115, 128, 137, 143, 144, 145, 146, 147, 149]

st = PDB.PDBParser(QUIET=True).get_structure("t", SD / "tnf_trimer.pdb")[0]
ShrakeRupley().compute(st, level="R")

# Numbering guard. These target files keep 1TNF's own numbering (6..157), so stored == mature.
# Only BoltzGen's OUTPUT renumbers to 1..152. An earlier version of this script assumed
# stored == mature-5 and silently reported every residue five positions off, which is the same
# off-by-five that gen_spec.py guards with a sequence assertion. Assert it here too.
from tnf_common import HUMAN as _H
_probe = {86: "S", 87: "Y", 88: "Q", 90: "K", 91: "V", 92: "N", 146: "E", 149: "Q"}
for _n, _want in _probe.items():
    _r = st["A"][_n]
    _got = THREE2ONE.get(_r.get_resname(), "X")
    assert _got == _want or _H[_n - 1] == _got, (
        f"numbering guard failed: A/{_n} is {_got}, expected {_want}")
print("numbering guard OK: stored residue numbers equal mature numbering\n")

heavy = np.array([a.coord for ch in st for r in ch if r.id[0] == " "
                  for a in r if a.element != "H"])

def res_at(ch: str, mature: int):
    """These files retain 1TNF numbering, so the mature number indexes directly."""
    try:
        return st[ch][mature]
    except KeyError:
        return None

info = {}
for ch in ("A", "B", "C"):
    for p in EPITOPE_ALL:
        r = res_at(ch, p)
        if r is None:
            continue
        atoms = [a.coord for a in r if a.element != "H"]
        if not atoms:
            continue
        cen = np.mean(atoms, axis=0)
        # outward normal: away from the whole-protein centroid
        normal = cen - heavy.mean(axis=0)
        n = normal / (np.linalg.norm(normal) + 1e-9)
        probe = cen + 5.0 * n           # a probe sitting just off the surface
        d = np.linalg.norm(heavy - probe, axis=1)
        concavity = int((d < 10.0).sum())    # atoms enclosing the probe
        aa = THREE2ONE.get(r.get_resname(), "X")
        info[(ch, p)] = {"cen": cen, "sasa": r.sasa, "aa": aa,
                         "concavity": concavity, "apolar": aa in APOLAR,
                         "divergent": p in DIVERGENT}

print("Per-residue designability of the TNFR2 epitope (trimer context)")
print(f"{'res':>8} {'aa':>3} {'SASA':>6} {'concavity':>10} {'apolar':>7} {'divergent':>10}")
for (ch, p), v in sorted(info.items(), key=lambda kv: -kv[1]["concavity"]):
    if ch != "A" or v["sasa"] < 10:
        continue
    print(f"{ch}/{p:<6} {v['aa']:>3} {v['sasa']:>6.1f} {v['concavity']:>10} "
          f"{str(v['apolar']):>7} {str(v['divergent']):>10}")

# Score tight clusters of 3-4 conserved exposed residues
print("\nBest TIGHT hotspot clusters (3-4 residues, conserved, exposed, <=12 A span):")
cands = [(k, v) for k, v in info.items()
         if not v["divergent"] and v["sasa"] >= 25 and k[0] in ("A", "B")]
import itertools
scored = []
for combo in itertools.combinations(cands, 3):
    pts = np.array([v["cen"] for _, v in combo])
    span = float(np.linalg.norm(pts[:, None] - pts[None], axis=-1).max())
    if span > 12.0:
        continue
    chains = {k[0] for k, _ in combo}
    conc = float(np.mean([v["concavity"] for _, v in combo]))
    apol = sum(1 for _, v in combo if v["apolar"])
    sasa = float(np.mean([v["sasa"] for _, v in combo]))
    # prefer concave, somewhat apolar, single-protomer, compact
    score = conc + 25 * apol - 8 * span - (40 if len(chains) > 1 else 0)
    scored.append((score, combo, span, conc, apol, sasa, chains))
scored.sort(reverse=True, key=lambda x: x[0])
seen = set()
shown = 0
for score, combo, span, conc, apol, sasa, chains in scored:
    key = frozenset(k for k, _ in combo)
    if key in seen:
        continue
    seen.add(key)
    label = ",".join(f"{k[0]}/{k[1]}{v['aa']}" for k, v in combo)
    print(f"  score={score:7.1f} span={span:4.1f}A concavity={conc:5.1f} apolar={apol} "
          f"meanSASA={sasa:5.1f} chains={sorted(chains)}  {label}")
    shown += 1
    if shown >= 10:
        break
