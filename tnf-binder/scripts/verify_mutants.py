"""Confirm each His mutant carries HIS at exactly the requested binder positions."""
import csv
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from make_his_mutants import binder_chain

variants, mut_dir, parents_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
rows = list(csv.DictReader(variants.open()))
by_name = {r["design"]: r for r in rows}

ok = bad = 0
for cif in sorted(mut_dir.glob("*.cif")):
    r = by_name.get(cif.stem)
    if r is None:
        continue
    want = [int(x) for x in r["his_positions"].split()]
    st = gemmi.read_structure(str(cif)); st.setup_entities()
    ch = binder_chain(st)
    seq = gemmi.one_letter_code([x.name for x in ch])
    # parent sequence for comparison
    psrc = next((p for p in parents_dir.rglob(f"{r['parent']}.cif")), None)
    pst = gemmi.read_structure(str(psrc)); pst.setup_entities()
    pseq = gemmi.one_letter_code([x.name for x in binder_chain(pst)])

    got_his = [i + 1 for i, c in enumerate(seq) if c == "H"]
    parent_his = [i + 1 for i, c in enumerate(pseq) if c == "H"]
    added = sorted(set(got_his) - set(parent_his))
    changed = [i + 1 for i, (a, b) in enumerate(zip(pseq, seq)) if a != b]
    good = (added == sorted(want)) and (changed == sorted(want)) and (r["sequence"] == seq)
    print(f"{cif.stem:28s} want={want} added_his={added} changed={changed} "
          f"seq_matches_csv={r['sequence'] == seq} {'OK' if good else 'MISMATCH'}")
    ok, bad = (ok + 1, bad) if good else (ok, bad + 1)

print(f"\n{ok} correct, {bad} mismatched")
