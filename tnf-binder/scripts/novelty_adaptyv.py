"""Reproduce Adaptyv/Proteinbase novelty scoring as closely as possible before submitting.

Their published scale (general proteins, 1-4; submissions must reach >=3):

  4  de novo        sequence identity <=30%  AND  less than moderate structural similarity
  3  partly novel   identity >30%  OR  moderate structural similarity -- but only one of them
  2  familiar       identity >70% alone, OR high structural similarity,
                    OR identity >30% together with moderate structural similarity
  1  known          identity >70% AND at least moderate structural similarity

  structural similarity is TM-score weighted by coverage:
     high      >=70% of the sequence covered by domains matching at TM-score >= 0.8
     moderate   >70% coverage at TM-score >= 0.5

Their pipeline uses MMseqs2 over SwissProt, PDB, patent sequences, a therapeutic-antibody
database and PLAbDab, plus ESMFold -> FoldSeek -> TM-align for the structural half.

This script covers the SEQUENCE half with MMseqs2 against the SwissProt and PDB-seqres copies
already on disk. It deliberately does NOT claim to reproduce the structural half; that needs
FoldSeek against a fold database. An exact-kmer screen (the earlier gate here) only catches
near-verbatim copies and would pass a design that still scores 2 on their scale.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

MMSEQS = Path.home() / "tools/mmseqs/bin/mmseqs"
# A hit must cover at least half the design to count as "sequence identity" in the sense the
# 30%/70% thresholds intend; short local alignments otherwise dominate a max-identity rule.
MIN_QUERY_COVERAGE = 0.5


def easy_search(query: Path, target: Path, workdir: Path) -> dict[str, tuple[float, str]]:
    """design -> (best percent identity, best hit id) against one database."""
    out = workdir / f"hits_{target.stem}.m8"
    cmd = [str(MMSEQS), "easy-search", str(query), str(target), str(out), str(workdir / "tmp"),
           "--format-output", "query,target,pident,alnlen,qlen,qcov,evalue,bits",
           "-s", "7.5", "--max-seqs", "300", "-e", "1e5", "--threads", "8"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not out.exists():
        print(f"  mmseqs failed on {target.name}: {r.stderr.strip()[-300:]}", file=sys.stderr)
        return {}
    # Rank by BITSCORE and require real coverage. Taking max percent-identity instead reports
    # tiny local alignments: it claimed 88.8% identity for a design whose genuine best hit is
    # 36% over 83% coverage, because a handful of residues matched well somewhere.
    best: dict[str, tuple[float, str]] = {}
    best_bits: dict[str, float] = {}
    for line in out.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        q, t = parts[0], parts[1]
        pid, qcov, bits = float(parts[2]), float(parts[5]), float(parts[7])
        pid = pid * 100 if pid <= 1.0 else pid        # mmseqs may report 0-1
        if qcov < MIN_QUERY_COVERAGE:
            continue
        if q not in best_bits or bits > best_bits[q]:
            best_bits[q] = bits
            best[q] = (pid, t)
    return best


def score(identity: float, structural: str) -> int:
    """structural is one of none|moderate|high (unknown treated as none, and flagged)."""
    high = structural == "high"
    moderate = structural in ("moderate", "high")
    if identity > 70 and moderate:
        return 1
    if identity > 70 or high or (identity > 30 and moderate):
        return 2
    if identity > 30 or moderate:
        return 3
    return 4


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", type=Path)
    ap.add_argument("--dbs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--structural", default="unknown",
                    choices=["unknown", "none", "moderate", "high"],
                    help="structural similarity verdict, if known from FoldSeek")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        per_db = {db.stem: easy_search(args.fasta, db, wd) for db in args.dbs}

    names = sorted({q for d in per_db.values() for q in d})
    rows = []
    for q in names:
        best_pid, best_hit, best_db = 0.0, "", ""
        for dbname, hits in per_db.items():
            if q in hits and hits[q][0] > best_pid:
                best_pid, best_hit, best_db = hits[q][0], hits[q][1], dbname
        struct = "none" if args.structural == "unknown" else args.structural
        s = score(best_pid, struct)
        rows.append({
            "design": q, "best_identity_pct": round(best_pid, 1),
            "best_hit": best_hit, "db": best_db,
            "structural_input": args.structural,
            "novelty_score_sequence_only": s,
            "passes_min_3": s >= 3,
        })

    print(f"{'design':34s} {'best_id%':>9s} {'score':>6s} {'pass>=3':>8s}  best hit")
    for r in sorted(rows, key=lambda x: -x["best_identity_pct"]):
        print(f"{r['design']:34s} {r['best_identity_pct']:9.1f} "
              f"{r['novelty_score_sequence_only']:6d} {str(r['passes_min_3']):>8s}  "
              f"{r['best_hit'][:44]}")

    if args.structural == "unknown":
        print("\nNOTE: structural similarity not measured, so these are an UPPER BOUND on the")
        print("      score. A design at <=30% identity still drops to 3 with moderate fold")
        print("      similarity, and to 2 with high fold similarity. FoldSeek is needed to settle it.")

    if args.out and rows:
        with args.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
