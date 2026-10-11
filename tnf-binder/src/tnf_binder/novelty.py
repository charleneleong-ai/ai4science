"""Search backends that feed the novelty verifiers.

Two independent questions, two tools, because they answer different things and only one of
them decided the real outcome:

  MMseqs2 over Swiss-Prot / PDB-seqres   sequence identity
  FoldSeek over CATH50                   fold similarity

Fold is the binding constraint. Two submitted designs differed by 0.6% sequence identity and
Proteinbase passed one and failed the other; the failure was a beta-sandwich matching the
immunoglobulin superfamily. A sequence screen cannot see that.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from tnf_binder.verifiers.fold_novelty import FoldHit

MIN_QUERY_COVERAGE = 0.5


@dataclass(frozen=True)
class SequenceSearch:
    """MMseqs2 over one or more FASTA corpora."""

    binary: Path
    databases: tuple[Path, ...]
    threads: int = 8

    def best_identity(self, fasta: Path) -> dict[str, tuple[float, str]]:
        """query -> (percent identity, best hit), ranked by bitscore.

        Ranked by bitscore with a coverage floor, never by maximum percent identity: that
        reports tiny local alignments and once claimed 88.8% identity for a design whose
        genuine best hit was 36% over 83% coverage.
        """
        best: dict[str, tuple[float, str]] = {}
        best_bits: dict[str, float] = {}
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for database in self.databases:
                hits = work / f"{database.stem}.m8"
                subprocess.run(
                    [str(self.binary), "easy-search", str(fasta), str(database), str(hits),
                     str(work / "tmp"), "--format-output",
                     "query,target,pident,alnlen,qlen,qcov,evalue,bits",
                     "-s", "7.5", "--max-seqs", "300", "-e", "1e5",
                     "--threads", str(self.threads)],
                    capture_output=True, text=True, check=False)
                if not hits.exists():
                    continue
                for line in hits.read_text().splitlines():
                    parts = line.split("\t")
                    if len(parts) < 8:
                        continue
                    query, target = parts[0], parts[1]
                    identity, coverage, bits = float(parts[2]), float(parts[5]), float(parts[7])
                    identity = identity * 100 if identity <= 1.0 else identity
                    if coverage < MIN_QUERY_COVERAGE:
                        continue
                    if query not in best_bits or bits > best_bits[query]:
                        best_bits[query], best[query] = bits, (identity, target)
        return best


@dataclass(frozen=True)
class FoldSearch:
    """FoldSeek over a fold database, returning hits for the fold-novelty verifier."""

    binary: Path
    database: Path
    threads: int = 8

    def hits(self, structure_dir: Path) -> dict[str, list[FoldHit]]:
        """design stem -> structural hits. Query with BINDER-ONLY structures.

        A whole complex would match TNF-alpha itself and report a fold hit that says nothing
        about the design.
        """
        found: dict[str, list[FoldHit]] = {}
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            out = work / "hits.m8"
            subprocess.run(
                [str(self.binary), "easy-search", str(structure_dir), str(self.database),
                 str(out), str(work / "tmp"), "--format-output",
                 "query,target,alntmscore,qcov,tcov,evalue,bits",
                 "--alignment-type", "1", "-e", "10", "--max-seqs", "500",
                 "--threads", str(self.threads)],
                capture_output=True, text=True, check=False)
            if not out.exists():
                return found
            for line in out.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) < 7:
                    continue
                stem = parts[0].replace("_binder.pdb", "").replace(".pdb", "")
                try:
                    found.setdefault(stem, []).append(
                        FoldHit(float(parts[2]), float(parts[3]), parts[1]))
                except ValueError:
                    continue
        return found


def kmer_corpus(flat_files: tuple[Path, ...]):
    """Exact k-mer lookup via grep -oF, for SequenceNoveltyVerifier.

    One pass per corpus using Aho-Corasick, rather than holding hundreds of millions of
    k-mers in memory.
    """
    def lookup(kmers: list[str]) -> set[str]:
        with tempfile.TemporaryDirectory() as tmp:
            patterns = Path(tmp) / "patterns.txt"
            patterns.write_text("\n".join(kmers) + "\n")
            hits: set[str] = set()
            for flat in flat_files:
                result = subprocess.run(["grep", "-oF", "-f", str(patterns), str(flat)],
                                        capture_output=True, text=True, check=False)
                hits |= {line for line in result.stdout.splitlines() if line}
            return hits
    return lookup
