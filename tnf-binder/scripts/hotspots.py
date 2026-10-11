"""Pick the compact binder epitope: max conserved SASA, min divergent SASA, must contain a conserved basic."""
import sys
from pathlib import Path
import numpy as np
from Bio import PDB

SASA = {20:68.4,21:168.6,23:165.3,31:192.9,32:91.4,33:86.0,63:3.6,67:120.1,71:94.6,72:71.9,
        73:47.8,75:45.8,77:34.6,85:78.4,86:80.7,87:102.0,88:154.5,89:72.2,90:108.4,91:42.3,
        92:56.4,97:31.4,111:34.9,113:42.3,115:56.8,128:144.3,137:70.0,143:14.5,144:83.1,
        145:87.6,146:53.1,147:55.1,149:0.0}
DIVERGENT = {20,31,71,72,73,85,89,97,111}
BASIC_CONSERVED = {32, 90, 128}
RADIUS = 11.0   # ~ footprint of a 60-120 aa minibinder

st = PDB.PDBParser(QUIET=True).get_structure("x", Path(sys.argv[1]) / "1tnf.pdb")[0]

def cen(ch, num):
    try:
        return np.mean([a.coord for a in st[ch][num] if a.element != "H"], axis=0)
    except KeyError:
        return None

# every epitope residue on chains A and B (the receptor cleft pair)
pts = [(ch, p, cen(ch, p)) for ch in ("A", "B") for p in SASA if cen(ch, p) is not None]

print(f"Compact epitope candidates (sphere r={RADIUS} A centred on each conserved exposed residue)\n")
rows = []
for ch0, p0, c0 in pts:
    if p0 in DIVERGENT or SASA[p0] < 40:
        continue
    near = [(ch, p) for ch, p, c in pts if np.linalg.norm(c - c0) <= RADIUS]
    cons_a = sum(SASA[p] for ch, p in near if p not in DIVERGENT)
    div_a  = sum(SASA[p] for ch, p in near if p in DIVERGENT)
    basics = sorted({p for ch, p in near if p in BASIC_CONSERVED})
    chains = {ch for ch, p in near}
    rows.append((cons_a - 2 * div_a, ch0, p0, cons_a, div_a, basics, near, chains))

rows.sort(reverse=True, key=lambda r: r[0])
for score, ch0, p0, ca, da, basics, near, chains in rows[:8]:
    frac = 100 * ca / (ca + da)
    print(f"centre {ch0}/{p0}: score={score:6.0f}  conserved={ca:5.0f} A^2  divergent={da:5.0f} A^2 "
          f"({frac:.0f}% conserved)")
    print(f"   pH anchors: {basics or 'NONE'}   protomers: {sorted(chains)}")
    print(f"   residues: {sorted(near, key=lambda x: (x[0], x[1]))}\n")
