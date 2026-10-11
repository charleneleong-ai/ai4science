"""Report which TNF-alpha residues a design actually contacts, handling merged target chains.

BindCraft concatenates the target protomers into one chain, so aligning that chain to a single
copy of human TNF maps only the first protomer and silently mis-scores everything after it.
This splits a merged chain into 157-residue protomer blocks before mapping to mature numbering.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import PDB

sys.path.insert(0, str(Path(__file__).parent))
from score_designs import load, classify_chains, chain_seq
from tnf_common import HUMAN, divergent_positions

EPITOPE_POCKET = {113, 115, 143, 144, 146, 147}
EPITOPE_LOBE = {86, 87, 88, 90, 91, 92, 33}
DIV = divergent_positions()
L = len(HUMAN)


def protomer_blocks(ch) -> dict:
    """target residue object -> mature resnum, splitting concatenated protomers."""
    res = [r for r in ch if r.id[0] == " "]
    seq = chain_seq(ch)
    out = {}
    # Find each protomer start by locating where the chain matches the start of HUMAN.
    starts = [i for i in range(len(seq) - 20)
              if seq[i:i + 20] == HUMAN[:20] or seq[i:i + 20] in HUMAN[:40]]
    if not starts:
        # fall back: assume contiguous copies of length L offset to 1TNF's first residue
        starts = list(range(0, len(res), L))
    # offset of the first modelled residue within HUMAN
    off = HUMAN.find(seq[:20]) if seq[:20] in HUMAN else 0
    for i, r in enumerate(res):
        blk = i // (len(res) // max(1, round(len(res) / L)) or len(res))
        within = i - blk * (len(res) // max(1, round(len(res) / L)) or len(res))
        out[r] = off + within + 1
    return out


for path in [Path(p) for p in sys.argv[1:]]:
    m = load(path)
    t, b = classify_chains(m)
    print(f"\n===== {path.name}")
    print(f"  targets={[c.id for c in t]} ({[len(chain_seq(c)) for c in t]} res) "
          f"binders={[c.id for c in b]} ({[len(chain_seq(c)) for c in b]} res)")
    if not t or not b:
        continue
    tnum = {}
    for ch in t:
        tnum.update(protomer_blocks(ch))
    tatoms = [a for ch in t for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    ns = PDB.NeighborSearch(tatoms)
    hits = Counter()
    for ch in b:
        for r in (x for x in ch if x.id[0] == " "):
            for a in (x for x in r if x.element != "H"):
                for ta in ns.search(a.coord, 4.5):
                    n = tnum.get(ta.get_parent())
                    if n:
                        hits[n] += 1
    tot = sum(hits.values())
    top = hits.most_common(20)
    print(f"  total contact atom-pairs={tot}")
    print(f"  top contacted mature residues: {[(n, c) for n, c in top]}")
    print(f"  of pocket epitope {sorted(EPITOPE_POCKET)}: hit {sorted(set(hits) & EPITOPE_POCKET)}")
    print(f"  of lobe epitope   {sorted(EPITOPE_LOBE)}: hit {sorted(set(hits) & EPITOPE_LOBE)}")
    div = sum(c for n, c in hits.items() if n in DIV)
    print(f"  divergent-residue contact fraction: {div / tot:.3f}" if tot else "  no contacts")
