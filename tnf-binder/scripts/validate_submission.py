"""Pre-flight check on submission CSVs against Proteinbase's stated requirements.

Requirements: columns name, sequence, molecule_class (first three); molecule_class in
{protein, nanobody, scfv, fab_kappa, fab_lambda}; length 10-250; de novo and unique.
Ranked with the top design first, so row order is meaningful.
"""
import csv
import sys
from pathlib import Path

VALID = {"protein", "nanobody", "scfv", "fab_kappa", "fab_lambda"}
STD = set("ACDEFGHIKLMNPQRSTVWY")
LIMITS = {"tnf_track1.csv": 40, "tnf_track2.csv": 20, "tnf_track3.csv": 20}

root = Path(sys.argv[1])
all_ok = True
for p in sorted(root.glob("tnf_track*.csv")):
    rows = list(csv.DictReader(p.open()))
    problems = []
    if not rows:
        problems.append("empty file")
    else:
        hdr = list(rows[0].keys())[:3]
        if hdr != ["name", "sequence", "molecule_class"]:
            problems.append(f"first three columns are {hdr}")
        names = [r["name"] for r in rows]
        seqs = [r["sequence"] for r in rows]
        if len(set(names)) != len(names):
            problems.append("duplicate names")
        if len(set(seqs)) != len(seqs):
            problems.append("duplicate sequences")
        cap = LIMITS.get(p.name)
        if cap and len(rows) > cap:
            problems.append(f"{len(rows)} rows exceeds cap {cap}")
        for r in rows:
            if r["molecule_class"] not in VALID:
                problems.append("bad molecule_class " + str(r["molecule_class"]))
                break
            bad = set(r["sequence"]) - STD
            if bad:
                problems.append("non-standard residues " + str(sorted(bad)) + " in " + r["name"])
                break
            if not (10 <= len(r["sequence"]) <= 250):
                problems.append("length " + str(len(r["sequence"])) + " out of range in " + r["name"])
                break
        lens = [len(s) for s in seqs] or [0]
        verdict = "OK" if not problems else "; ".join(problems)
        print(p.name + ": " + str(len(rows)) + " rows, lengths " +
              str(min(lens)) + "-" + str(max(lens)) + " -> " + verdict)
    all_ok = all_ok and not problems

print("ALL FILES VALID" if all_ok else "VALIDATION FAILED")
sys.exit(0 if all_ok else 1)
