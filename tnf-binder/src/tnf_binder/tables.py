"""The CSVs each pipeline stage hands to the next, and the two coercions they all needed.

Every stage writes a CSV and the next stage reads it, so the same two rules were re-stated in
every script: a missing file is an empty table rather than an error, because a stage that has
not run yet is the normal case; and a blank or malformed numeric cell falls back to a default
rather than raising, because BoltzGen and PROPKA both leave cells empty when a metric could
not be computed, and one such cell must not sink a whole batch.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def read_rows(path: Path | None) -> list[dict[str, str]]:
    """Rows of a CSV, or nothing at all if the stage that writes it has not run."""
    if path is None or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle))


def number(row: Mapping[str, Any], key: str, default: float = 0.0) -> float:
    """A numeric cell, with blanks and junk coerced to `default` instead of raising."""
    try:
        return float(row.get(key) or default)
    except (TypeError, ValueError):
        return default


def write_rows(path: Path, rows: Sequence[Mapping[str, Any]],
               columns: Sequence[str] | None = None) -> None:
    """Write rows to a CSV, keyed on `columns` if given, else on the first row's keys."""
    fieldnames = list(columns) if columns else list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
