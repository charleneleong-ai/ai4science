"""Is this design novel, or did the generator reproduce something that already exists?

`precedent.py` asks whether a coordination motif has been *seen before* and treats a match as
reassuring. This asks the opposite question about the whole chain, and treats a match as
disqualifying: a generator that emits a natural protein has produced nothing, however well it
scores. The top-scoring construct of one campaign — best interface metric, lowest predicted
error, an ideally-placed histidine — was verbatim human ubiquitin, and no geometry tier can
see that, because its geometry is genuinely excellent.

Two independent questions, two tools:

    SequenceSearch   MMseqs2 over sequence corpora (Swiss-Prot, PDB seqres, ...)
    FoldSearch       FoldSeek over a structure database (CATH50, PDB, ...)

**Fold is the binding constraint, not sequence.** Two designs submitted to the same reviewer
differed by 0.6 percentage points of sequence identity and one passed while the other failed;
the failure was an immunoglobulin-like beta-sandwich at TM 0.893. A sequence screen cannot see
that, and an exact k-mer screen is blind to it by construction.

Both wrappers shell out and return empty on failure rather than raising, so a missing binary
degrades to "cannot judge" at the caller rather than taking down a batch.
"""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

# Rank sequence hits by bitscore above this query coverage, never by maximum percent identity.
MIN_QUERY_COVERAGE = 0.5


@dataclass(frozen=True)
class FoldHit:
    """One structural hit: TM-score, the fraction of the query it covers, and what it matched."""

    tm_score: float
    coverage: float
    target: str


@dataclass(frozen=True)
class SequenceHit:
    """One sequence hit, with the coverage that makes its identity meaningful."""

    identity_pct: float
    coverage: float
    bits: float
    target: str


@dataclass(frozen=True)
class SequenceSearch:
    """MMseqs2 over one or more FASTA corpora."""

    binary: Path
    databases: tuple[Path, ...]
    threads: int = 8
    min_coverage: float = MIN_QUERY_COVERAGE

    def best(self, fasta: Path) -> dict[str, SequenceHit]:
        """query id -> its best hit, ranked by **bitscore** above the coverage floor.

        Never by maximum percent identity. That reports whichever tiny local alignment happens
        to match best and once claimed 88.8% identity for a design whose genuine best hit was
        36% over 83% coverage — which would have scored a de novo design as memorised.
        """
        best: dict[str, SequenceHit] = {}
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for database in self.databases:
                hits = work / f"{database.stem}.m8"
                subprocess.run(
                    [
                        str(self.binary),
                        "easy-search",
                        str(fasta),
                        str(database),
                        str(hits),
                        str(work / "tmp"),
                        "--format-output",
                        "query,target,pident,alnlen,qlen,qcov,evalue,bits",
                        "-s",
                        "7.5",
                        "--max-seqs",
                        "300",
                        "-e",
                        "1e5",
                        "--threads",
                        str(self.threads),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if hits.exists():
                    self.absorb(hits.read_text(), best)
        return best

    def absorb(self, tabular: str, best: dict[str, SequenceHit]) -> None:
        """Merge one MMseqs2 tabular block into `best`, keeping the highest bitscore per query."""
        for line in tabular.splitlines():
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            try:
                identity, coverage, bits = (
                    float(parts[2]),
                    float(parts[5]),
                    float(parts[7]),
                )
            except ValueError:
                continue
            # MMseqs reports identity as a fraction or a percentage depending on version
            identity = identity * 100 if identity <= 1.0 else identity
            if coverage < self.min_coverage:
                continue
            query = parts[0]
            if query not in best or bits > best[query].bits:
                best[query] = SequenceHit(identity, coverage, bits, parts[1])


@dataclass(frozen=True)
class FoldSearch:
    """FoldSeek over a structure database."""

    binary: Path
    database: Path
    threads: int = 8

    def hits(self, structure_dir: Path) -> dict[str, list[FoldHit]]:
        """design stem -> structural hits. Query with **binder-only** structures.

        A whole complex matches its own target and reports a fold hit that says nothing about
        the design.
        """
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            out = work / "hits.m8"
            subprocess.run(
                [
                    str(self.binary),
                    "easy-search",
                    str(structure_dir),
                    str(self.database),
                    str(out),
                    str(work / "tmp"),
                    "--format-output",
                    "query,target,alntmscore,qcov,tcov,evalue,bits",
                    "--alignment-type",
                    "1",
                    "-e",
                    "10",
                    "--max-seqs",
                    "500",
                    "--threads",
                    str(self.threads),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            return parse_foldseek(out.read_text()) if out.exists() else {}


def parse_foldseek(tabular: str) -> dict[str, list[FoldHit]]:
    """FoldSeek tabular output -> hits per design stem."""
    found: dict[str, list[FoldHit]] = {}
    for line in tabular.splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        stem = (
            parts[0].replace("_binder.pdb", "").replace(".pdb", "").replace(".cif", "")
        )
        try:
            hit = FoldHit(float(parts[2]), float(parts[3]), parts[1])
        except ValueError:
            continue
        found.setdefault(stem, []).append(hit)
    return found


def kmer_corpus(flat_files: Sequence[Path]) -> Callable[[Sequence[str]], set[str]]:
    """An exact k-mer lookup over flat sequence files, via one `grep -oF` pass per corpus.

    grep's `-f` is Aho-Corasick: all patterns are matched in a single pass, so this stays
    linear in corpus size and never holds hundreds of millions of k-mers in memory. Exact
    matching is the point — it is decisive where an alignment score is arguable, and a design
    sharing a 20-mer with a natural protein is not de novo whatever its identity percentage.
    """

    def lookup(kmers: Sequence[str]) -> set[str]:
        if not kmers:
            return set()
        with tempfile.TemporaryDirectory() as tmp:
            patterns = Path(tmp) / "patterns.txt"
            patterns.write_text("\n".join(kmers) + "\n")
            hits: set[str] = set()
            for flat in flat_files:
                result = subprocess.run(
                    ["grep", "-oF", "-f", str(patterns), str(flat)],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                hits |= {line for line in result.stdout.splitlines() if line}
            return hits

    return lookup
