"""Map the human TNF-alpha epitope onto mouse numbering for BindCraft2 dual-target design.

Objective 2 wants one binder for both orthologs. BindCraft2 can optimise both interfaces at
once, but each target needs hotspots in ITS OWN residue numbering. The human construct is
157 aa and the mouse 156 aa with a gap in the alignment, so positions do not simply transfer --
this derives them, and asserts the mapped residues are the identical ones we rely on.
"""
import sys
from pathlib import Path

from Bio import Align

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN, MOUSE

# conserved pocket (what we target) and the lobe, in human mature numbering
POCKET = [113, 115, 143, 144, 146, 147]
LOBE = [86, 87, 88, 90, 91, 92]
# human/mouse divergent positions to forbid contact with (BindCraft2 "coldspots")
DIVERGENT = [20, 31, 71, 72, 73, 85, 89, 97, 111]

al = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)
ha, ma = (str(x) for x in al.align(HUMAN, MOUSE)[0])

h2m, hi, mi = {}, 0, 0
for hc, mc in zip(ha, ma):
    if hc != "-":
        hi += 1
    if mc != "-":
        mi += 1
    if hc != "-" and mc != "-":
        h2m[hi] = (mi, hc, mc)

print(f"{'human':>6} {'aa':>3} -> {'mouse':>6} {'aa':>3}  identical?")
mouse_pocket, mouse_lobe = [], []
for p in POCKET + LOBE:
    if p not in h2m:
        print(f"{p:6d}  --   (no aligned mouse position)")
        continue
    m, hc, mc = h2m[p]
    same = "yes" if hc == mc else "NO"
    print(f"{p:6d} {hc:>3} -> {m:6d} {mc:>3}  {same}")
    assert hc == mc, f"human {p}{hc} maps to mouse {m}{mc}: not conserved, should not be a hotspot"
    (mouse_pocket if p in POCKET else mouse_lobe).append(m)

mouse_cold = [h2m[p][0] for p in DIVERGENT if p in h2m]

print(f"\nhuman pocket hotspots : {POCKET}")
print(f"mouse pocket hotspots : {mouse_pocket}")
print(f"human lobe hotspots   : {LOBE}")
print(f"mouse lobe hotspots   : {mouse_lobe}")
print(f"human coldspots (divergent): {DIVERGENT}")
print(f"mouse coldspots (divergent): {mouse_cold}")

# The human target PDB keeps 1TNF numbering (6..157); the mouse model was predicted from the
# 156-aa soluble sequence and is numbered from 1.
print("\nNOTE: human PDB uses 1TNF numbering (mature == stored); the mouse prediction starts at 1,")
print("      so mouse numbers above are already the model's own numbering.")
