"""Build design-ready TNF target files and verify the A/91 hotspot patch geometry."""
import sys
from pathlib import Path
import numpy as np
from Bio import PDB
from Bio.PDB import Select

SD = Path(sys.argv[1])
HOTSPOT = [("A", p) for p in (86, 87, 88, 90, 91, 92)] + [("B", 33), ("B", 147)]

st = PDB.PDBParser(QUIET=True).get_structure("tnf", SD / "1tnf.pdb")[0]

# --- completeness check: 1TNF is a 1989 2.6 A structure, expect gaps ---
print("Residue coverage per chain in 1TNF (mature numbering 1-157):")
present = {}
for ch in st:
    nums = sorted(r.id[1] for r in ch if r.id[0] == " ")
    present[ch.id] = set(nums)
    missing = [n for n in range(1, 158) if n not in present[ch.id]]
    print(f"  chain {ch.id}: {len(nums)} res, range {nums[0]}-{nums[-1]}, missing {len(missing)}: {missing}")

print("\nHotspot residues present in all chains?")
for ch, p in HOTSPOT:
    print(f"  {ch}/{p}: {'OK' if p in present.get(ch, ()) else 'MISSING'}")

def cen(ch, num):
    try:
        return np.mean([a.coord for a in st[ch][num] if a.element != "H"], axis=0)
    except KeyError:
        return None

hs = np.array([c for ch, p in HOTSPOT if (c := cen(ch, p)) is not None])
print(f"\nHotspot patch: n={len(hs)}  max internal span={np.linalg.norm(hs[:,None]-hs[None],axis=-1).max():.1f} A")
print(f"  centroid={hs.mean(0).round(1)}")

# is protomer C far enough to truncate away?
cmin = min(np.linalg.norm(a.coord - hs, axis=-1).min() for a in st["C"].get_atoms() if a.element != "H")
print(f"\nMin distance chain C -> hotspot patch: {cmin:.1f} A "
      f"({'C is remote, A+B truncation is SAFE' if cmin > 12 else 'C is close, KEEP FULL TRIMER'})")

# divergent residues near the patch, for the design-time exclusion list
DIV = {20:"P/H",31:"R/Q",71:"S/D",72:"T/-",73:"H/Y",85:"V/I",89:"T/E",97:"I/V",111:"A/L"}
print("\nDivergent residues within 10 A of hotspot patch (avoid heavy contact):")
for p, lab in DIV.items():
    ds = [np.linalg.norm(c - hs, axis=-1).min() for ch in ("A","B") if (c := cen(ch, p)) is not None]
    if ds and min(ds) < 10:
        print(f"  {p} ({lab}): {min(ds):.1f} A")

class Keep(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, c): return c.id in self.chains
    def accept_residue(self, r): return r.id[0] == " "

io = PDB.PDBIO(); io.set_structure(st)
for name, chains in (("tnf_trimer.pdb", "ABC"), ("tnf_AB.pdb", "AB")):
    io.save(str(SD / name), Keep(chains))
    print(f"wrote {name}")
