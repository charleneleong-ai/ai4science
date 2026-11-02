"""Fold novelty: the check the organisers actually apply, which sequence identity cannot do.

`novelty.py` rejects a design that is the target's own chain handed back, by sequence. That is
necessary and not sufficient. Two facts make it insufficient here:

- The generator's recorded failure on this protocol is collapse onto a natural protein, not onto
  the target: the TNF campaign saw 2 of 8 designs come back as verbatim ubiquitin, which scores
  near zero identity against every chain of a pMHC and sails through a target-only gate.
- The competition's own novelty check is fold-led. A beta-sandwich design of ours failed it at
  37.5% sequence identity while a helical design passed at 38.1%, which is recorded in
  `get_cath50.sh` on the design host. Beta-2 microglobulin is an immunoglobulin beta-sandwich,
  so the designs most at risk are the ones that are beta-2 microglobulin shaped at a sequence
  identity low enough for `novelty.py` to clear them.

So this runs FoldSeek over the designed chain alone, against CATH50 at domain level and against
the PDB, and reports the best structural match. The binder chain has to be extracted first: a
search with the whole co-fold pose matches the HLA trivially and says nothing about the design.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import gemmi
import typer
from common import binder_chain, load, target_sequences, write_rows

app = typer.Typer(add_completion=False)

# FoldSeek's alignment TM-score. 0.5 is the conventional line for "same fold", but on this
# cohort it does not separate anything: all 11 designs of the first batch matched something at
# TM 0.71-0.93, because a small helical bundle resembles one of the thousands already deposited.
# A TM-only reject would therefore throw away every design. What does separate them is the
# sequence identity to the structural match: 95.1% to 5xs3 chain B is a known protein handed
# back, while 0.0-16.2% at the same TM is a de novo bundle that happens to share a common fold.
SAME_FOLD = 0.50
COPIED_IDENTITY = 0.30

# CATH's immunoglobulin-like beta-sandwich, which is beta-2 microglobulin's own architecture and
# the class a submission of ours failed the organisers' novelty check in, at 37.5% sequence
# identity, while a helical design passed at 38.1% (recorded in get_cath50.sh on the design
# host). Only the AlphaFold-derived CATH50 entries carry the superfamily code in their name, so
# this catches the hits it can read and the rest are left to the identity test above.
IG_SANDWICH = "2.60.40"

FORMAT = "query,target,fident,alntmscore,evalue,alnlen"


@dataclass(frozen=True)
class FoldHit:
    design: str
    database: str
    target: str
    identity: float
    tm_score: float
    evalue: float
    length: int

    @property
    def verdict(self) -> str:
        """Inadmissible, admissible-but-unremarkable, or structurally distinct.

        Deliberately not a TM-only verdict: see the note on SAME_FOLD. `copy_of_known` is the
        only value that is a reject on its own.
        """
        if self.tm_score >= SAME_FOLD and self.identity >= COPIED_IDENTITY:
            return "copy_of_known"
        if IG_SANDWICH in self.target:
            return "ig_sandwich"
        if self.tm_score >= SAME_FOLD:
            return "shared_fold"
        return "distinct_fold"

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "database": self.database,
            "closest_fold": self.target,
            "sequence_identity": round(self.identity, 4),
            "tm_score": round(self.tm_score, 4),
            "evalue": self.evalue,
            "alignment_length": self.length,
            "verdict": self.verdict,
        }


def extract_binder(structure: gemmi.Structure, out: Path, chain: str) -> Path:
    extracted = gemmi.Structure()
    extracted.add_model(gemmi.Model("1"))
    copied = gemmi.Chain("A")
    for residue in structure[0][chain]:
        copied.add_residue(residue)
    extracted[0].add_chain(copied)
    extracted.setup_entities()
    out.parent.mkdir(parents=True, exist_ok=True)
    extracted.write_pdb(str(out))
    return out


def run_foldseek(binaries: Path, query_dir: Path, database: Path, out: Path, tmp: Path) -> Path:
    """One easy-search over the whole directory; FoldSeek is far faster batched than per file."""
    tmp.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(binaries / "foldseek"), "easy-search",
        str(query_dir), str(database), str(out), str(tmp),
        "--format-output", FORMAT,
        "--alignment-type", "1",
        "-e", "10000",
        "--max-seqs", "1000",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise typer.BadParameter(
            f"foldseek failed ({result.returncode}): {result.stderr.strip()[-800:]}"
        )
    return out


def best_hits(table: Path, database: str) -> dict[str, FoldHit]:
    """Highest TM-score per design; FoldSeek sorts by e-value, which is not what gates here.

    The highest-TM hit is what defines a design's fold, and picking anything else was measured
    to be wrong. Ranking instead by how disqualifying a hit is -- so that a lower-TM
    immunoglobulin match outranked a strong helical one -- labelled 8 of 11 designs
    immunoglobulin beta-sandwich at TM 0.42-0.72, including one whose own best matches were
    TM 0.865 and 0.995 to something else. This search runs at `-e 10000` over 1,000 hits per
    query on purpose, so some 2.60.40 entry appears for nearly every design; reading the tail
    of that list as a verdict reports the most spurious hit rather than the most informative.
    """
    best: dict[str, FoldHit] = {}
    with table.open() as handle:
        for line in handle:
            if not (parts := line.rstrip("\n").split("\t")):
                continue
            if len(parts) < 6:
                continue
            query, target, fident, tm, evalue, alnlen = parts[:6]
            design = Path(query).stem
            hit = FoldHit(
                design=design,
                database=database,
                target=target,
                identity=float(fident),
                tm_score=float(tm),
                evalue=float(evalue),
                length=int(alnlen),
            )
            if design not in best or hit.tm_score > best[design].tm_score:
                best[design] = hit
    return best


@app.command()
def main(
    poses: Path = typer.Argument(..., help="directory of co-fold or generator poses (.cif/.pdb)"),
    out: Path = typer.Option(Path("fold_novelty.csv")),
    work: Path = typer.Option(Path("fold_work"), help="extracted binder chains and FoldSeek tmp"),
    binaries: Path = typer.Option(Path.home() / "tools/foldseek/bin"),
    databases: str = typer.Option(
        f"cath50={Path.home()}/tools/cathdb/cath50,pdb={Path.home()}/tools/fsdb/pdb"
    ),
    references: Path = typer.Option(
        Path("targets/novelty_references.json"), help="target sequences, to exclude the target"
    ),
) -> None:
    exclude = target_sequences(references)
    found = sorted(p for p in poses.rglob("*.cif") if p.is_file())
    if not found:
        raise typer.BadParameter(f"{poses} holds no .cif poses")

    chains_dir = work / "binders"
    for pose in found:
        structure = load(pose)
        chain = binder_chain(pose, structure, exclude)
        extract_binder(structure, chains_dir / f"{pose.stem}.pdb", chain)
    typer.echo(f"extracted {len(found)} binder chains into {chains_dir}")

    rows: list[FoldHit] = []
    for entry in databases.split(","):
        name, _, path = entry.partition("=")
        table = run_foldseek(binaries, chains_dir, Path(path), work / f"{name}.m8", work / "tmp")
        hits = best_hits(table, name)
        typer.echo(f"\n=== {name}: {len(hits)} of {len(found)} designs matched anything ===")
        for hit in sorted(hits.values(), key=lambda h: -h.tm_score)[:10]:
            typer.echo(
                f"  {hit.verdict:10s} {hit.design:24s} TM={hit.tm_score:.3f}"
                f" id={100 * hit.identity:5.1f}%  {hit.target}"
            )
        counts: dict[str, int] = {}
        for hit in hits.values():
            counts[hit.verdict] = counts.get(hit.verdict, 0) + 1
        unmatched = len(found) - len(hits)
        typer.echo(f"  {counts}, plus {unmatched} with no hit at all")
        rows.extend(hits.values())

    written = write_rows(rows, out, label="hits")
    typer.echo(f"\nwrote {out} ({written} rows)")


if __name__ == "__main__":
    app()
