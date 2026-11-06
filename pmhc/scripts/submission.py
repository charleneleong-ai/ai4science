"""Build and validate the Challenge 5 submission from the measurements, not by hand.

Proteinbase's format, as the live submit page states it -- the page and the FAQ disagreed and
the page won, a lesson this repo paid for on an earlier track and which
`tnf_binder.submission` records in full:

  molecule_class  one of single_chain / nanobody / scfv / fab_kappa / fab_lambda. "protein" is
                  NOT accepted; these de novo miniproteins are single_chain.
  cap             20 sequences per submission.
  length          10-250 aa per chain.

The rules are **duplicated** from `tnf_binder.submission`, deliberately and under deadline.
tnf-binder is not a dependency of this project, but it does sit in this same repo and does
import cleanly, so the honest statement is that importing it was not attempted rather than that
it was impossible: it would mean a path dependency from one campaign into another's source
tree, when these are portal rules that belong to neither. The real fix is a shared module both
campaigns depend on -- `glp1r` is a third track against the same portal -- and it is recorded in
common.py's follow-up note so that the next format change has somewhere to look. What this copy
loses is tnf-binder's regression anchor: a test that re-validates the CSV actually uploaded.

**Ranking.** By mean peptide-dependence drop over repeated diffusion samples, never by a single
draw. The pilot screen ranked on one sample per condition and put two designs in its top four
whose drops reverse sign when measured five times -- which is the whole reason this reads a
multi-sample table. `drop_significant` says whether a design's drop clears its own sampling
noise and `versus_control` whether it is distinguishable from the one protein known to do the
job; both come from `cofold_followup.comparison`, which owns that rule because it is the module
that sees `n`. Neither is a filter: the cap is 20, the field is smaller than that after the
gates, and a submission of 20 with 3 above noise is reported as exactly that.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import typer
from cofold_followup import SIGMA, comparison
from cofold_specs import CONTROL, design_sequences
from common import rounded, rows_by, write_rows, wrote

app = typer.Typer(add_completion=False)

MOLECULE_CLASS = "single_chain"
VALID_CLASSES = frozenset({MOLECULE_CLASS, "nanobody", "scfv", "fab_kappa", "fab_lambda"})
SUBMISSION_COLUMNS = ["name", "sequence", "molecule_class"]
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH, MAX_LENGTH = 10, 250
MAX_ROWS = 20

PREFIX = "pmhc"


def row_problem(row: Mapping[str, str]) -> str | None:
    """The first way one submission row violates Proteinbase's stated format, if any."""
    sequence = row["sequence"]
    if row["molecule_class"] not in VALID_CLASSES:
        return f"bad class {row['molecule_class']}"
    if set(sequence) - STANDARD_AA:
        return f"non-standard residues in {row['name']}"
    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        return f"length {len(sequence)} in {row['name']}"
    return None


def csv_problems(rows: Sequence[Mapping[str, str]]) -> list[str]:
    """Every format violation in one submission file.

    A wrong header returns alone: `row_problem` indexes the columns it checks, so running it
    over rows that do not have them raises KeyError from inside the validator instead of
    reporting the column problem already found. `check` reads files this code did not write.
    """
    if not rows:
        return ["empty"]
    if (columns := list(rows[0])) != SUBMISSION_COLUMNS:
        return [f"columns {columns}, expected {SUBMISSION_COLUMNS}"]
    problems = []
    if len({r["sequence"] for r in rows}) != len(rows):
        problems.append("duplicate sequences")
    if len({r["name"] for r in rows}) != len(rows):
        problems.append("duplicate names")
    if len(rows) > MAX_ROWS:
        problems.append(f"{len(rows)} rows > {MAX_ROWS} limit")
    return problems + [p for r in rows if (p := row_problem(r))]


def problems_in(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    """Read a submission CSV back off disk and check it. Used to verify what was just written
    as well as to check a file from elsewhere -- the point is to validate the bytes that get
    uploaded, including the header `write_rows` emitted and any quoting the writer applied."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    return rows, csv_problems(rows)


@dataclass(frozen=True)
class Control:
    """The 7OW6 receptor's own drop -- the same for every candidate, so it gets one name."""

    drop: float
    sem: float

    @classmethod
    def read(cls, row: Mapping[str, str]) -> Control:
        return cls(drop=float(row["drop_mean"]), sem=float(row["drop_sem"]))


@dataclass(frozen=True)
class SubmissionRow:
    """Exactly the three columns the portal accepts, and nothing else."""

    name: str
    sequence: str

    @property
    def row(self) -> dict[str, object]:
        return {"name": self.name, "sequence": self.sequence, "molecule_class": MOLECULE_CLASS}


@dataclass(frozen=True)
class Candidate:
    rank: int
    design: str
    sequence: str
    drop_mean: float
    drop_sd: float
    drop_sem: float
    sem_basis: str
    drop_own_sd: float | None
    control: Control
    evidence: Mapping[str, str]

    @property
    def name(self) -> str:
        return f"{PREFIX}_{self.rank:02d}"

    @property
    def drop_significant(self) -> bool:
        """Whether this design's own drop clears the standard error of its own mean."""
        return comparison(self.drop_mean, self.drop_sem, 0.0, 0.0) == "above"

    @property
    def versus_control(self) -> str:
        """Above, tied with, or below the one molecule known to discriminate this epitope."""
        return comparison(self.drop_mean, self.drop_sem, self.control.drop, self.control.sem)

    @property
    def submission(self) -> SubmissionRow:
        return SubmissionRow(name=self.name, sequence=self.sequence)

    @property
    def row(self) -> dict[str, object]:
        return {
            "rank": self.rank,
            "name": self.name,
            "source_design": self.design,
            "variant": self.evidence.get("variant", ""),
            "length": len(self.sequence),
            "drop_mean": round(self.drop_mean, 4),
            "drop_sd": round(self.drop_sd, 4),
            "drop_sem": round(self.drop_sem, 4),
            "sem_basis": self.sem_basis,
            "drop_own_sd": rounded(self.drop_own_sd),
            "drop_significant": self.drop_significant,
            "versus_control": self.versus_control,
            "peptide_share": self.evidence.get("peptide_share", ""),
            "mutation_contacts": self.evidence.get("mutation_contacts", ""),
            "copied_fraction": self.evidence.get("copied_fraction", ""),
            "fold_verdict": self.evidence.get("fold_verdict", ""),
        }


def candidates(
    spread: Mapping[str, Mapping[str, str]],
    sequences: Mapping[str, str],
    manifest: Mapping[str, Mapping[str, str]],
    count: int,
) -> list[Candidate]:
    if (control_row := spread.get(CONTROL)) is None:
        raise typer.BadParameter(
            f"no {CONTROL} row in the spread table: the drops would have no scale and must not"
            " be ranked. Re-run the co-fold pass with the control's specs included."
        )
    control = Control.read(control_row)

    ranked = sorted(
        (row for design, row in spread.items() if design != CONTROL),
        key=lambda row: -float(row["drop_mean"]),
    )[:count]
    if not ranked:
        raise typer.BadParameter("the spread table holds the control and no designs")

    # Every missing sequence at once, named in one message -- the shape cofold_specs.py uses for
    # the same invariant, which it establishes before any spec is written.
    if absent := sorted(r["design"] for r in ranked if r["design"] not in sequences):
        raise typer.BadParameter(f"no sequence for selected design(s) {absent}")

    return [
        Candidate(
            rank=rank,
            design=row["design"],
            sequence=sequences[row["design"]],
            drop_mean=float(row["drop_mean"]),
            drop_sd=float(row["drop_sd"]),
            drop_sem=float(row["drop_sem"]),
            sem_basis=row.get("sem_basis", ""),
            drop_own_sd=(
                float(own) if (own := row.get("drop_own_sd")) not in (None, "") else None
            ),
            control=control,
            evidence=manifest.get(row["design"], {}),
        )
        for rank, row in enumerate(ranked, start=1)
    ]


@app.command()
def build(
    spread: Path = typer.Option(
        Path("cofold_spread_pooled.csv"),
        help="written by cofold_followup pool-replicates, or report-spread for a single pass"
    ),
    manifest: Path = typer.Option(Path("cofold_manifest.csv")),
    gen_out: Path = typer.Option(Path("gen_out")),
    out: Path = typer.Option(Path("submission/pmhc_submission.csv")),
    evidence: Path = typer.Option(Path("submission/pmhc_submission_evidence.csv")),
    count: int = typer.Option(MAX_ROWS),
) -> None:
    # The spread table first: its control row is the cheapest and most fundamental check, and
    # parsing every generator metrics CSV before making it is work thrown away on a bad run.
    rows = rows_by(spread, "design")
    if CONTROL not in rows:
        raise typer.BadParameter(f"no {CONTROL} row in {spread}; the drops have no scale")
    chosen = candidates(rows, design_sequences(gen_out), rows_by(manifest, "design"), count)

    typer.echo(
        f"  {'name':10s} {'design':20s} {'len':>4s} {'drop':>8s} {'sem':>6s}"
        f" {'>noise':>7s} {'vs ctrl':>8s}"
    )
    for c in chosen:
        typer.echo(
            f"  {c.name:10s} {c.design:20s} {len(c.sequence):4d}"
            f" {c.drop_mean:+8.3f} {c.drop_sem:6.3f}"
            f" {'yes' if c.drop_significant else 'no':>7s} {c.versus_control:>8s}"
        )

    verdicts = [c.versus_control for c in chosen]
    tally = {v: verdicts.count(v) for v in ("above", "tied", "below")}
    control = chosen[0].control
    typer.echo(
        f"\n{len(chosen)} submitted; {sum(c.drop_significant for c in chosen)} have a drop"
        f" clearing the standard error of their own mean at {SIGMA:g} sigma"
    )
    typer.echo(
        f"against the control's {control.drop:+.3f} +- {control.sem:.3f}: {tally['above']}"
        f" above, {tally['tied']} tied, {tally['below']} below"
    )

    typer.echo(wrote(out, write_rows([c.submission for c in chosen], out)))
    typer.echo(wrote(evidence, write_rows(chosen, evidence)))

    _, problems = problems_in(out)
    if problems:
        typer.echo("FORMAT PROBLEMS: " + "; ".join(problems))
        raise typer.Exit(1)
    typer.echo("format OK: columns, classes, residues, lengths, uniqueness, cap")


@app.command()
def check(path: Path = typer.Argument(..., help="a submission CSV to validate")) -> None:
    rows, problems = problems_in(path)
    typer.echo(f"{path}: {len(rows)} rows")
    for problem in problems:
        typer.echo(f"  PROBLEM {problem}")
    if problems:
        raise typer.Exit(1)
    lengths = sorted(len(r["sequence"]) for r in rows)
    typer.echo(f"format OK; lengths {lengths[0]}-{lengths[-1]} aa")


if __name__ == "__main__":
    app()
