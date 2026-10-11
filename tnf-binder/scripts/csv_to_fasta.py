"""Dump a CSV column pair to FASTA, for the novelty gate."""
import csv
import sys
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])
name_col = sys.argv[3] if len(sys.argv) > 3 else "spec"
seq_col = sys.argv[4] if len(sys.argv) > 4 else "sequence"

rows = list(csv.DictReader(src.open()))
with out.open("w") as f:
    n = 0
    for r in rows:
        s = (r.get(seq_col) or "").strip()
        if s:
            f.write(f">{r[name_col]}\n{s}\n")
            n += 1
print(f"wrote {n} sequences -> {out}")
