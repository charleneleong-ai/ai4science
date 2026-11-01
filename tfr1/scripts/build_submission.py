"""Turn the selection table into the two files an entry consists of, and refuse a bad one.

The upload wants `name,sequence,molecule_class` and nothing else, so everything that justifies a
row -- its switch, which chassis it came from, what the independent co-fold saw -- goes in a
companion evidence table that travels with the methodology instead.

The format rules are imported from `tnf_binder.submission` rather than restated. They were read
off the live submit page and then verified against three accepted uploads, including the one that
established that `molecule_class` is `single_chain` where the FAQ says `protein`; a second copy
of a rule that subtle is a copy that will be wrong. That package is a sibling in this repo and
not installed here, so its source directory is put on the path explicitly.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

import typer
from csvio import interface_reads, number, read_rows, write_rows

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tnf-binder" / "src"))

from tnf_binder.submission import (  # noqa: E402  (needs the path insert above)
    MAX_LENGTH,
    MIN_LENGTH,
    MOLECULE_CLASS,
    SUBMISSION_COLUMNS,
    csv_problems,
    write_submission,
)

# Columns carried through unchanged from the selection, and the ones lifted from the independent
# co-fold under a `chai_` prefix. Named once so the row and the header cannot disagree.
CARRIED = ("parent", "epitope", "substitutions", "switch_log_units", "fold_weaker_at_ph6",
           "coupled_sites", "epitope_conservation", "novelty_score", "novelty_note")
CO_FOLD = ("interface_residues", "hotspots_touched", "nearest_hotspot_angstroms",
           "pose_agreement")
EVIDENCE_COLUMNS = ["rank", "name", "source_design", "length", *CARRIED, "parent_iptm",
                    *(f"chai_{field}" for field in CO_FOLD)]

app = typer.Typer(add_completion=False)


def problems(rows: list[dict[str, str]], candidates: list[dict[str, str]]) -> list[str]:
    """Every format violation, plus the one check the shared rules cannot make.

    `csv_problems` knows the upload's own rules -- the columns, the cap, the legal class, the
    length range, duplicate names and sequences. What it cannot see is whether a row still
    matches the table it was built from, so the only check added here is that cross-reference:
    a sequence silently truncated or re-wrapped between selection and here would pass every
    format rule and submit a different protein than the evidence table describes.
    """
    return csv_problems(rows) + [
        f"{row['name']}: sequence is {len(row['sequence'])} aa, table says {candidate['length']}"
        for row, candidate in zip(rows, candidates, strict=True)
        if len(row["sequence"]) != int(candidate["length"])]


@app.command()
def main(
    candidates: Annotated[Path, typer.Option(help="the selected designs, in rank order")],
    out: Annotated[Path, typer.Option(help="the CSV to upload")],
    evidence: Annotated[Path, typer.Option(help="what justifies each row, for the methodology")],
    interfaces: Annotated[Path | None, typer.Option(
        help="pose_interface output, to carry the independent co-fold read per design")] = None,
    prefix: Annotated[str, typer.Option(help="submission name stem")] = "tfr1",
) -> None:
    chosen = read_rows(candidates)
    # The independent co-fold specifically, never the generator's own pose: both are in that
    # table and only one of them belongs under a `chai_` column.
    seen = interface_reads(interfaces, "chai")

    rows, evidence_rows = [], []
    for rank, candidate in enumerate(chosen, start=1):
        name = f"{prefix}_{rank:02d}"
        sequence = candidate["sequence"].strip().upper()
        rows.append({"name": name, "sequence": sequence, "molecule_class": MOLECULE_CLASS})
        read = seen.get(candidate["name"], {})
        evidence_rows.append(
            {"rank": rank, "name": name, "source_design": candidate["name"],
             "length": len(sequence),
             "parent_iptm": round(number(candidate["parent_iptm"]), 3)}
            | {field: candidate[field] for field in CARRIED}
            | {f"chai_{field}": read.get(field, "") for field in CO_FOLD})

    if found := problems(rows, chosen):
        for problem in found:
            typer.echo(f"REJECTED  {problem}", err=True)
        raise typer.Exit(1)

    write_submission(out, SUBMISSION_COLUMNS,
                     [[row[column] for column in SUBMISSION_COLUMNS] for row in rows])
    write_rows(evidence, evidence_rows, EVIDENCE_COLUMNS)
    lengths = [len(row["sequence"]) for row in rows]
    typer.echo(f"{len(rows)} rows, {min(lengths)}-{max(lengths)} aa "
               f"(legal range {MIN_LENGTH}-{MAX_LENGTH}), class {MOLECULE_CLASS}")
    typer.echo(f"wrote {out} and {evidence}")
    if not seen:
        typer.echo("no independent co-fold read carried; the evidence columns for it are empty")


if __name__ == "__main__":
    app()
