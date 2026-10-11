"""Fold-similarity half of the Proteinbase novelty score, via FoldSeek + TM-align.

This exists because sequence novelty is NOT the binding constraint. A submitted design at 37.5%
sequence identity FAILED their novelty check while another at 38.1% passed; the failing one was
a beta-sandwich (FN3-like), the passing one a helical bundle. Their pipeline runs ESMFold ->
FoldSeek -> TM-align, and scores structural similarity as TM-score weighted by coverage:

  high      >=70% of the sequence covered by domains matching at TM-score >= 0.8
  moderate   >70% coverage at TM-score >= 0.5

Combined with sequence identity, their 1-4 scale is:
  4 de novo   identity <=30% and less than moderate structural similarity
  3           identity >30% XOR moderate structural similarity
  2           identity >70%, or high structural similarity, or (identity >30% and moderate)
  1           identity >70% and at least moderate structural similarity

Searched against CATH50 (domain-level, matching their domain-segmentation step) rather than the
full PDB set, because the volume this runs on is at 99%.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

FOLDSEEK = Path.home() / "tools/foldseek/bin/foldseek"


def search(query_dir: Path, db: Path, workdir: Path) -> dict[str, list[tuple[float, float, str]]]:
    """design -> [(tm-score, query coverage, target), ...]"""
    out = workdir / "fold_hits.m8"
    cmd = [str(FOLDSEEK), "easy-search", str(query_dir), str(db), str(out), str(workdir / "tmp"),
           "--format-output", "query,target,alntmscore,qcov,tcov,evalue,bits",
           "--alignment-type", "1", "-e", "10", "--max-seqs", "500", "--threads", "8"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not out.exists():
        print(f"foldseek failed: {r.stderr.strip()[-400:]}", file=sys.stderr)
        return {}
    hits: dict[str, list[tuple[float, float, str]]] = {}
    for line in out.read_text().splitlines():
        p = line.split("\t")
        if len(p) < 7:
            continue
        q, t = p[0], p[1]
        try:
            tm, qcov = float(p[2]), float(p[3])
        except ValueError:
            continue
        hits.setdefault(q, []).append((tm, qcov, t))
    return hits


def structural_verdict(hits: list[tuple[float, float, str]]) -> tuple[str, float, float, str]:
    """(none|moderate|high, best tm, its coverage, target) per their coverage-weighted rule."""
    best_high = max(((tm, c, t) for tm, c, t in hits if tm >= 0.8 and c >= 0.70),
                    default=None, key=lambda x: x[0])
    if best_high:
        return "high", best_high[0], best_high[1], best_high[2]
    best_mod = max(((tm, c, t) for tm, c, t in hits if tm >= 0.5 and c > 0.70),
                   default=None, key=lambda x: x[0])
    if best_mod:
        return "moderate", best_mod[0], best_mod[1], best_mod[2]
    best_any = max(hits, default=(0.0, 0.0, ""), key=lambda x: x[0])
    return "none", best_any[0], best_any[1], best_any[2]


def novelty_score(identity: float, structural: str) -> int:
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
    ap.add_argument("binder_dir", type=Path, help="dir of binder-only PDBs")
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--identities", type=Path,
                    help="CSV from novelty_adaptyv.py: design,best_identity_pct")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    ident: dict[str, float] = {}
    if args.identities and args.identities.exists():
        for r in csv.DictReader(args.identities.open()):
            ident[r["design"]] = float(r["best_identity_pct"])

    with tempfile.TemporaryDirectory() as td:
        hits = search(args.binder_dir, args.db, Path(td))

    rows = []
    for q in sorted(hits):
        verdict, tm, cov, tgt = structural_verdict(hits[q])
        stem = q.replace("_binder.pdb", "").replace(".pdb", "")
        pid = next((v for k, v in ident.items() if k in stem or stem in k), float("nan"))
        s = novelty_score(pid if pid == pid else 0.0, verdict)
        rows.append({"design": stem, "identity_pct": round(pid, 1) if pid == pid else "",
                     "structural": verdict, "best_tm": round(tm, 3),
                     "tm_coverage": round(cov, 3), "closest_fold": tgt,
                     "novelty_score": s, "passes_min_3": s >= 3})

    print(f"{'design':44s} {'id%':>6s} {'fold':>9s} {'TM':>6s} {'cov':>6s} "
          f"{'score':>6s} {'pass':>6s}  closest")
    for r in rows:
        print(f"{r['design']:44s} {str(r['identity_pct']):>6s} {r['structural']:>9s} "
              f"{r['best_tm']:6.3f} {r['tm_coverage']:6.3f} {r['novelty_score']:6d} "
              f"{str(r['passes_min_3']):>6s}  {r['closest_fold'][:26]}")

    if args.out and rows:
        with args.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
