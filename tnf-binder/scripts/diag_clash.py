"""Is the clash count real interpenetration, or an artefact of how the file is written?"""
import sys
from pathlib import Path

import numpy as np
from Bio import PDB

sys.path.insert(0, str(Path(__file__).parent))
from score_designs import load, classify_chains, chain_seq

for p in [Path(x) for x in sys.argv[1:]]:
    m = load(p)
    t, b = classify_chains(m)
    print(f"\n===== {p.parent.name}/{p.name}")
    for ch in m:
        res = [r for r in ch if r.id[0] == " "]
        atoms = [a for r in res for a in r]
        role = "TARGET" if ch in t else "BINDER"
        print(f"   {role} ch={ch.id} res={len(res)} atoms={len(atoms)} "
              f"atoms/res={len(atoms)/max(1,len(res)):.1f}")
    ta = [a for ch in t for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    ba = [a for ch in b for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    if not ta or not ba:
        continue
    T = np.array([a.coord for a in ta]); B = np.array([a.coord for a in ba])
    # pairwise min distances
    d = np.linalg.norm(B[:, None, :] - T[None, :, :], axis=-1)
    print(f"   global min interatomic dist = {d.min():.2f} A")
    for cut in (1.0, 1.5, 2.0, 2.2, 2.5, 3.0):
        print(f"     pairs < {cut:.1f} A: {(d < cut).sum()}")
    # which atom names are involved in the sub-2.2 contacts?
    i, j = np.where(d < 2.2)
    if len(i):
        names = {}
        for x, y in list(zip(i, j))[:400]:
            names[(ba[x].get_id(), ta[y].get_id())] = names.get((ba[x].get_id(), ta[y].get_id()), 0) + 1
        top = sorted(names.items(), key=lambda kv: -kv[1])[:8]
        print(f"   most common clashing atom-name pairs: {top}")
