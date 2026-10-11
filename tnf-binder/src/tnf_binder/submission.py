"""Proteinbase's submission format, as the live submit page states it.

The page and the FAQ disagreed, and the page won:

  molecule_class  one of single_chain / nanobody / scfv / fab_kappa / fab_lambda. "protein" is
                  NOT accepted; these de novo miniproteins are single_chain.
  cap             20 sequences per SUBMISSION, not the 40/20/20 per-track split the FAQ
                  described.
  length          10-250 aa per chain.

Novelty is auto-scored after upload and must reach 3/4, so a memorised sequence is actively
harmful rather than merely a wasted slot -- which is why the novelty verifiers gate before this
and not after.

Owned here rather than in the CLI because the assembly scripts write the files the CLI then
validates; two statements of the same rules would let a file pass the writer and fail the
checker.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from tnf_binder.design import STANDARD_AA

VALID_CLASSES = frozenset({"single_chain", "nanobody", "scfv", "fab_kappa", "fab_lambda"})
MOLECULE_CLASS = "single_chain"
SUBMISSION_COLUMNS = ["name", "sequence", "molecule_class"]
MIN_LENGTH, MAX_LENGTH = 10, 250
MAX_ROWS = 20


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


def csv_problems(rows: list[dict[str, str]], max_rows: int = MAX_ROWS) -> list[str]:
    """Every format violation in one submission file."""
    if not rows:
        return ["empty"]
    problems = []
    columns = list(rows[0])[:3]
    if columns != SUBMISSION_COLUMNS:
        problems.append(f"columns {columns}")
    if len({r["sequence"] for r in rows}) != len(rows):
        problems.append("duplicate sequences")
    if len({r["name"] for r in rows}) != len(rows):
        problems.append("duplicate names")
    if len(rows) > max_rows:
        problems.append(f"{len(rows)} rows > {max_rows} limit")
    return problems + [p for r in rows if (p := row_problem(r))]


def write_submission(path: Path, columns: Sequence[str],
                     rows: Sequence[Sequence[Any]]) -> None:
    """Write a submission CSV with LF line endings, matching the official template exactly.

    csv.writer defaults to CRLF (RFC 4180); the template is LF-only and their parser has not
    been tested here, so match the template rather than the standard.
    """
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(list(columns))
        writer.writerows(rows)
