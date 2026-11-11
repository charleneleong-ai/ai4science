"""Reject designs that are the target's own chains handed back.

The first two designs this pipeline produced were beta-2 microglobulin, a chain of the target
complex, at better than 90% identity. The generator had been given the HLA heavy chain with its
beta-2 microglobulin interface vacant, and it filled the vacancy with the protein that belongs
there. Those are high-scoring poses of a known complex, not binders, and the challenge requires
designs to be de novo and zero-shot and forbids modifying an existing binder, so submitting one
risks disqualification rather than a low score.

A regurgitated chain can come back truncated or with scattered substitutions, so the gate is the
fraction of the design that is a verbatim copy of some target chain, measured by local alignment.
It is a hard reject, not a ranking term.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import typer
from Bio.Align import PairwiseAligner
from common import target_references, write_rows, wrote

app = typer.Typer(add_completion=False)

DEFAULT_THRESHOLD = 0.30


def local_aligner() -> PairwiseAligner:
    """Smith-Waterman, which is what sees a chain handed back in pieces.

    This gate used to score a global alignment's matches over the shorter sequence, and that
    cannot represent a truncated copy at all: measured on verbatim slices of the deposition, a
    40-residue piece of the HLA heavy chain scored 0.000 and a 60-residue piece of the TCR beta
    chain 0.000, both passed as novel, against a threat model that names truncation explicitly.
    The same slices score 1.000 here. Mismatches have to cost something or a local alignment
    extends through them for free.
    """
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 1.0
    a.mismatch_score = -1.0
    a.open_gap_score = -11.0
    a.extend_gap_score = -1.0
    return a


def global_aligner() -> PairwiseAligner:
    """The superseded whole-chain measure, still reported so the change stays auditable."""
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 1.0
    a.mismatch_score = 0.0
    a.open_gap_score = -1.0
    a.extend_gap_score = -0.5
    return a


def copied_fraction(query: str, reference: str, align: PairwiseAligner) -> float:
    """Identities in the best local alignment over the design's own length.

    Normalised by the design, not by the shorter sequence: against the 10-residue epitope a
    reference-normalised score would call four chance matches 40% copied.
    """
    if not query or not reference:
        return 0.0
    alignments = align.align(query, reference)
    if not len(alignments):
        return 0.0
    return alignments[0].counts().identities / len(query)


def global_identity(query: str, reference: str, align: PairwiseAligner) -> float:
    if not query or not reference:
        return 0.0
    return max(0.0, float(align.score(query, reference)) / min(len(query), len(reference)))


@dataclass(frozen=True)
class Verdict:
    design: str
    length: int
    worst_match: str
    copied: float
    whole_chain: float
    admissible: bool

    @property
    def row(self) -> dict[str, str | int | float | bool]:
        return {
            "design": self.design,
            "length": self.length,
            "closest_target_chain": self.worst_match,
            "copied_fraction": round(self.copied, 4),
            "whole_chain_identity": round(self.whole_chain, 4),
            "novel": self.admissible,
        }


def judge(
    designs: dict[str, str], refs: dict[str, str], threshold: float
) -> list[Verdict]:
    local, whole = local_aligner(), global_aligner()
    verdicts: list[Verdict] = []
    for name, sequence in designs.items():
        worst, key = max(
            (copied_fraction(sequence, ref, local), ref_name)
            for ref_name, ref in refs.items()
        )
        verdicts.append(
            Verdict(
                design=name,
                length=len(sequence),
                worst_match=key,
                copied=worst,
                whole_chain=max(
                    global_identity(sequence, ref, whole) for ref in refs.values()
                ),
                admissible=worst < threshold,
            )
        )
    return verdicts


def read_designs(metrics: Path, name_column: str, sequence_column: str) -> dict[str, str]:
    with metrics.open() as handle:
        return {
            row[name_column]: row[sequence_column]
            for row in csv.DictReader(handle)
            if row.get(sequence_column)
        }


@app.command()
def main(
    metrics: Path = typer.Argument(..., help="a BoltzGen *_metrics*.csv"),
    targets: Path = typer.Option(Path("targets/novelty_references.json")),
    out: Path = typer.Option(Path("novelty.csv")),
    threshold: float = typer.Option(
        DEFAULT_THRESHOLD, help="reject at or above this copied fraction"
    ),
    name_column: str = typer.Option("id"),
    sequence_column: str = typer.Option("designed_sequence"),
) -> None:
    designs = read_designs(metrics, name_column, sequence_column)
    if not designs:
        raise typer.BadParameter(f"{metrics} has no usable {sequence_column} values")
    verdicts = judge(designs, target_references(targets), threshold)

    for v in sorted(verdicts, key=lambda v: -v.copied):
        mark = "novel " if v.admissible else "REJECT"
        typer.echo(
            f"  {mark} {v.design:20s} len={v.length:3d}"
            f" {100 * v.copied:5.1f}% copied from {v.worst_match}"
            f" (whole-chain {100 * v.whole_chain:5.1f}%)"
        )
    rejected = sum(1 for v in verdicts if not v.admissible)
    typer.echo(
        f"{len(verdicts) - rejected} of {len(verdicts)} are novel at the {threshold:.0%} threshold"
    )

    typer.echo(wrote(out, write_rows(verdicts, out)))


if __name__ == "__main__":
    app()
