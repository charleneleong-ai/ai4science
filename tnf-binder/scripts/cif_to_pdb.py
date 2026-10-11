"""Convert mmCIF complexes to PDB, because propka reads PDB only."""
import sys
from pathlib import Path

import gemmi

out_dir = Path(sys.argv[1])
out_dir.mkdir(parents=True, exist_ok=True)
n = 0
for arg in sys.argv[2:]:
    for src in (sorted(Path(arg).rglob("*.cif")) if Path(arg).is_dir() else [Path(arg)]):
        st = gemmi.read_structure(str(src))
        st.setup_entities()
        st.remove_alternative_conformations()
        dest = out_dir / (src.stem + ".pdb")
        st.write_pdb(str(dest))
        n += 1
print(f"converted {n} structures -> {out_dir}")
