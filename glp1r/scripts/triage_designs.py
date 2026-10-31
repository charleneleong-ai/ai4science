"""Sequence-level triage, applied before anything expensive is spent on a design.

Every filter here is one the TNF-alpha campaign paid to learn, and each rejects
a specific observed failure rather than expressing a preference:

- **Composition degeneracy.** Poly-alanine and near-poly-glycine designs appeared
  at the top of a shortlist ranked on interface confidence. They are artefacts of
  the generator, not weak candidates, so they are rejected rather than ranked.
- **Free cysteine.** 18 of 20 designs in the TNF submission carried a single
  unpaired cysteine at the same position, inherited from one chassis, and it was
  only noticed at packaging time. A lone thiol risks disulfide scrambling,
  covalent dimerisation and aggregation.
- **Expression.** On 2,214 wet-lab-tested public designs the best predictor of
  expression found was glutamate fraction, within-target AUROC 0.724 -- better
  than any structural score tested, including ipTM at 0.607. It costs nothing.

Thresholds are declared here and are not tuned against results.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

# Competition limits.
MIN_LENGTH, MAX_LENGTH = 10, 250

# Composition gate, from the TNF campaign's observed failures.
MAX_SINGLE_RESIDUE_FRACTION = 0.30
MIN_DISTINCT_RESIDUES = 12
MAX_ALA_GLY_FRACTION = 0.45

# Expression: higher glutamate fraction predicted expression at AUROC 0.724.
# This is a *rank* signal, not a hard gate -- it is reported, not used to reject,
# because a threshold on it was never validated.
GLU_REFERENCE = 0.283  # the helix_0 chassis, which the predictor rates favourably


@dataclass(frozen=True)
class Verdict:
    name: str
    sequence: str
    problems: tuple[str, ...]
    glu_fraction: float
    cysteines: int

    @property
    def passes(self) -> bool:
        return not self.problems


def problems_with(sequence: str) -> tuple[str, ...]:
    found: list[str] = []
    length = len(sequence)
    if not sequence:
        return ("empty",)
    if not MIN_LENGTH <= length <= MAX_LENGTH:
        found.append(f"length {length} outside {MIN_LENGTH}-{MAX_LENGTH}")
    if non_standard := sorted(set(sequence) - STANDARD):
        found.append(f"non-standard residues {''.join(non_standard)}")

    counts = Counter(sequence)
    residue, top = counts.most_common(1)[0]
    if top / length > MAX_SINGLE_RESIDUE_FRACTION:
        found.append(f"{residue} is {top / length:.0%} of the sequence")
    if len(counts) < MIN_DISTINCT_RESIDUES:
        found.append(f"only {len(counts)} distinct residues")
    if (counts["A"] + counts["G"]) / length > MAX_ALA_GLY_FRACTION:
        found.append(f"Ala+Gly {(counts['A'] + counts['G']) / length:.0%}")
    return tuple(found)


def judge(name: str, sequence: str) -> Verdict:
    sequence = sequence.strip().upper()
    return Verdict(
        name=name,
        sequence=sequence,
        problems=problems_with(sequence),
        glu_fraction=sequence.count("E") / len(sequence) if sequence else 0.0,
        cysteines=sequence.count("C"),
    )


def read_sequences(path: Path) -> list[tuple[str, str]]:
    """Accept a FASTA or a CSV carrying name and sequence columns."""
    if path.suffix.lower() in {".fasta", ".fa"}:
        out: list[tuple[str, str]] = []
        name, chunks = None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(chunks)))
                name, chunks = line[1:].strip().split()[0], []
            elif line.strip():
                chunks.append(line.strip())
        if name:
            out.append((name, "".join(chunks)))
        return out
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    key = next(k for k in rows[0] if k.lower() in {"sequence", "seq"})
    label = next((k for k in rows[0] if k.lower() in {"name", "id", "design"}), None)
    return [(r[label] if label else f"design_{i:03d}", r[key]) for i, r in enumerate(rows)]


app = typer.Typer(add_completion=False)


@app.command()
def main(
    sequences: Annotated[Path, typer.Argument(help="FASTA or CSV of designs")],
    out: Annotated[Path, typer.Option(help="per-design verdict CSV")],
) -> None:
    verdicts = [judge(name, seq) for name, seq in read_sequences(sequences)]
    kept = [v for v in verdicts if v.passes]

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "length", "passes", "problems", "glu_fraction",
                         "cysteines", "sequence"])
        for v in verdicts:
            writer.writerow([v.name, len(v.sequence), "yes" if v.passes else "no",
                             "; ".join(v.problems), f"{v.glu_fraction:.3f}",
                             v.cysteines, v.sequence])

    reasons = Counter(p.split()[0] for v in verdicts for p in v.problems)
    typer.echo(f"{len(verdicts)} designs -> {len(kept)} pass ({len(kept) / max(len(verdicts), 1):.0%})")
    if reasons:
        typer.echo("rejections by first word of reason:")
        for reason, n in reasons.most_common():
            typer.echo(f"  {n:4d}  {reason}")
    if kept:
        free_cys = sum(1 for v in kept if v.cysteines == 1)
        glu = sorted(v.glu_fraction for v in kept)
        typer.echo(f"\nof those passing: {free_cys} carry exactly one cysteine (unpaired thiol risk)")
        typer.echo(f"glutamate fraction median {glu[len(glu) // 2]:.3f} (helix_0 reference {GLU_REFERENCE:.3f})")
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
