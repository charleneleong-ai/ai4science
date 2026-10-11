"""Write binder-only PDBs for fold comparison.

FoldSeek must see the binder alone: querying with the whole complex would match TNF-alpha
itself and report a fold hit that says nothing about the design.
"""
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN


def is_tnf(seq: str) -> bool:
    return any(seq[i:i + 25] in HUMAN for i in range(0, max(1, len(seq) - 25), 5))


out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
n = 0
for arg in sys.argv[2:]:
    src = Path(arg)
    paths = sorted(src.rglob("*.cif")) + sorted(src.rglob("*.pdb")) if src.is_dir() else [src]
    for p in paths:
        st = gemmi.read_structure(str(p))
        st.setup_entities()
        st.remove_alternative_conformations()
        keep = []
        for ch in st[0]:
            seq = gemmi.one_letter_code([r.name for r in ch])
            if seq and len(seq) < 260 and not is_tnf(seq):
                keep.append(ch.name)
        if not keep:
            continue
        for model in st:
            for name in [c.name for c in model]:
                if name not in keep:
                    model.remove_chain(name)
        st.setup_entities()
        dest = out / f"{p.stem}_binder.pdb"
        st.write_pdb(str(dest))
        n += 1
print(f"wrote {n} binder-only structures -> {out}")
