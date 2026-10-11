"""Find the largest spatially-contiguous, species-conserved, exposed patch in the TNFR2 epitope."""
import sys
from pathlib import Path
import numpy as np
from Bio import PDB

CONSERVED = [21,23,32,33,63,67,75,77,86,87,88,90,91,92,113,115,128,137,143,144,145,146,147,149]
DIVERGENT = [20,31,71,72,73,85,89,97,111]
EXPOSED = {21:168.6,23:165.3,32:91.4,33:86.0,63:3.6,67:120.1,75:45.8,77:34.6,86:80.7,87:102.0,
           88:154.5,90:108.4,91:42.3,92:56.4,113:42.3,115:56.8,128:144.3,137:70.0,143:14.5,
           144:83.1,145:87.6,146:53.1,147:55.1,149:0.0}

st = PDB.PDBParser(QUIET=True).get_structure("x", Path(sys.argv[1]) / "1tnf.pdb")[0]

def centroid(ch: str, num: int) -> np.ndarray | None:
    try:
        return np.mean([a.coord for a in st[ch][num] if a.element != "H"], axis=0)
    except KeyError:
        return None

# TNFR2 binds at the cleft between two adjacent protomers. Build the composite
# site from chain A plus its neighbour B, keeping only exposed conserved residues.
site: list[tuple[str, int, np.ndarray]] = []
for ch in ("A", "B"):
    for p in CONSERVED:
        if EXPOSED.get(p, 0) < 25:
            continue
        c = centroid(ch, p)
        if c is not None:
            site.append((ch, p, c))

coords = np.array([c for _, _, c in site])
labels = [f"{ch}/{p}" for ch, p, _ in site]
D = np.linalg.norm(coords[:, None] - coords[None], axis=-1)

# single-linkage cluster at 12 A to find contiguous surface patches
n = len(site)
parent = list(range(n))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i in range(n):
    for j in range(i + 1, n):
        if D[i, j] < 12.0:
            a, b = find(i), find(j)
            if a != b: parent[a] = b

clusters: dict[int, list[int]] = {}
for i in range(n):
    clusters.setdefault(find(i), []).append(i)

print("Contiguous conserved+exposed patches across the A/B protomer cleft:\n")
for k, idxs in sorted(clusters.items(), key=lambda kv: -len(kv[1])):
    if len(idxs) < 3:
        continue
    mem = [labels[i] for i in idxs]
    area = sum(EXPOSED[site[i][1]] for i in idxs)
    span = D[np.ix_(idxs, idxs)].max()
    chains = {site[i][0] for i in idxs}
    bridges = "BRIDGES BOTH PROTOMERS" if len(chains) > 1 else f"single protomer {chains}"
    print(f"patch n={len(idxs)}  SASA~{area:.0f} A^2  span={span:.1f} A  {bridges}")
    print(f"   {sorted(mem)}\n")

# distance from divergent residues to the best patch
best = max(clusters.values(), key=len)
bc = coords[best]
print("Min distance from each DIVERGENT (human/mouse differing) residue to best patch:")
for p in DIVERGENT:
    ds = [np.linalg.norm(centroid(ch, p) - bc, axis=-1).min()
          for ch in ("A", "B") if centroid(ch, p) is not None]
    if ds:
        print(f"  {p:>4}: {min(ds):>5.1f} A  {'<-- INSIDE patch, risk' if min(ds) < 8 else ''}")

print("\npH-switch anchors (conserved, exposed, basic) inside best patch:")
for ch, p, c in site:
    if p in (32, 90, 128) and any(np.allclose(c, coords[i]) for i in best):
        print(f"  {ch}/{p}")
