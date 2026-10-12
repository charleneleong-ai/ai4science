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

# Gates for a candidate competing for a free slot, fixed BEFORE the designs they judge were
# folded, so what survives is a pre-registered count and not a threshold fitted to the answer.
# Each is taken from a design already submitted rather than chosen to admit a particular number.
MIN_DELTA_Q = 0.47          # the weakest submitted design whose His is robust on every sample
MIN_MOUSE_RETAINED = 0.70   # the campaign's own floor; 9 of 12 submitted clear it, median 86%
MIN_ROBUST_HIS = 2          # at least two histidines hold the window on every diffusion sample
MAX_CLASHES = 0
# NOT pre-registered: the structural `flag` check below was added after reading the first
# geometry table, which showed two candidates at epitope_focus 0.00 -- bound somewhere else on
# the target entirely, and admitted by a clash count of zero. Declared here rather than folded
# in silently. It is not a threshold fitted to the result: `scripts/submit.py shortlist` has
# rejected on this same flag since the first campaign, and the omission was in failing to wire
# it into the newer gate.


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


def candidate_problem(record: Mapping[str, Any], *, max_clashes: int = MAX_CLASHES) -> str | None:
    """Why a measured candidate does not earn a slot, or None if it clears every gate.

    Separate from `row_problem`, which asks only whether a row is a legal upload. This asks
    whether the design is worth one of twenty slots, and it drops rather than down-ranks: the
    objectives are ordered pH-selectivity, then mouse cross-reactivity, then affinity, and a
    design that fails either of the first two is not rescued by being good at the third.

    Structural validity is checked first: a clashing or off-epitope complex is not a structure
    whose pKa and interface mean anything, so reporting it as a near miss on delta-Q would be
    reporting a measurement taken on the wrong object.
    """
    flag = str(record.get("flag") or "").strip().strip('"')
    # A clash flag is judged by the count, which the caller can calibrate; the other flags
    # (off_epitope, tiny_interface) are categorical and always disqualifying.
    other = [f for f in flag.split(",") if f and not f.startswith("clash(")]
    if other:
        return f"structural flag: {','.join(other)}"
    if float(record["clashes"] or 0) > max_clashes:
        return f"{record['clashes']} clashes"
    if float(record["delta_q_floor"]) < MIN_DELTA_Q:
        return f"delta_q {record['delta_q_floor']} < {MIN_DELTA_Q}"
    if float(record["mouse_retained"]) < MIN_MOUSE_RETAINED:
        return (f"mouse retained {float(record['mouse_retained']):.0%} "
                f"< {MIN_MOUSE_RETAINED:.0%}")
    robust = int(str(record["his_robust"]).split("/")[0])
    if robust < MIN_ROBUST_HIS:
        return f"only {record['his_robust']} His robust across samples"
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
