"""Flag designs that are memorised natural proteins rather than de novo.

The competition hard-filters on de novo design, and BoltzGen demonstrably collapses onto
natural folds -- the first batch returned verbatim human ubiquitin. A genuine de novo sequence
essentially never shares a long exact k-mer with a natural protein, so an exact 20-mer hit
against PDB seqres or Swiss-Prot is decisive evidence of memorisation.

Exact k-mer matching is used rather than a BLAST web service: it needs no API key or email
hand-off, has no rate limit, and for this yes/no question it is strictly more interpretable.
It will not catch a remote homologue -- it is a memorisation screen, not a homology search.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path


def read_sequences(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    if path.suffix == ".csv":
        for row in csv.DictReader(path.open()):
            s = (row.get("sequence") or "").strip().upper()
            if s:
                seqs[row.get("design") or row.get("name") or f"row{len(seqs)}"] = s
    else:  # fasta
        name, buf = None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if name and buf:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line.strip():
                buf.append(line.strip())
        if name and buf:
            seqs[name] = "".join(buf).upper()
    return seqs


def kmers(seq: str, k: int, stride: int) -> list[str]:
    if len(seq) < k:
        return []
    out = [seq[i:i + k] for i in range(0, len(seq) - k + 1, stride)]
    tail = seq[-k:]
    if tail not in out:
        out.append(tail)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("sequences", type=Path, help="FASTA or CSV with a 'sequence' column")
    ap.add_argument("--db", type=Path, required=True, help="dir holding *.flat reference files")
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    seqs = read_sequences(args.sequences)
    if not seqs:
        print("no sequences found")
        return

    kmer_owner: dict[str, set[str]] = {}
    for name, s in seqs.items():
        for km in kmers(s, args.k, args.stride):
            kmer_owner.setdefault(km, set()).add(name)
    print(f"{len(seqs)} sequences -> {len(kmer_owner)} distinct {args.k}-mers", file=sys.stderr)

    pat = args.db / "_patterns.txt"
    pat.write_text("\n".join(kmer_owner) + "\n")

    hits: set[str] = set()
    for flat in sorted(args.db.glob("*.flat")):
        # grep -oF with many patterns uses Aho-Corasick, so one pass per reference file.
        r = subprocess.run(["grep", "-oF", "-f", str(pat), str(flat)],
                           capture_output=True, text=True)
        found = {ln for ln in r.stdout.splitlines() if ln}
        print(f"  {flat.name}: {len(found)} matching k-mers", file=sys.stderr)
        hits |= found
    pat.unlink(missing_ok=True)

    rows = []
    for name, s in sorted(seqs.items()):
        matched = [km for km in kmers(s, args.k, args.stride) if km in hits]
        rows.append({
            "design": name,
            "length": len(s),
            "novel": not matched,
            "n_kmer_hits": len(matched),
            "example_hit": matched[0] if matched else "",
            "sequence": s,
        })

    n_novel = sum(1 for r in rows if r["novel"])
    print(f"\nnovel: {n_novel}/{len(rows)}   memorised: {len(rows) - n_novel}/{len(rows)}")
    for r in rows:
        if not r["novel"]:
            print(f"  REJECT {r['design']}: {r['n_kmer_hits']} hits, e.g. {r['example_hit']}")

    if args.out:
        with args.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
