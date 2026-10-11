"""Copy only the structurally sound designs forward, for His engineering.

his_engineer must never see junk: on a clashing poly-alanine backbone every position looks
"buried", and in testing it proposed 40-77 bogus His sites per design that way.
"""
import csv
import shutil
import sys
from pathlib import Path

scored, struct, dest = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
dest.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader(scored.open()))
keep = [r["design"] for r in rows if not (r.get("flag") or "").strip()]
n = 0
for d in keep:
    for cand in struct.rglob(f"{d}.*"):
        if cand.suffix in (".cif", ".pdb"):
            shutil.copy2(cand, dest / cand.name)
            n += 1
            break
print(f"   {n}/{len(rows)} designs passed scoring -> {dest}")
