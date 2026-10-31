"""The design table's reader and writer, in one place.

Three scripts in this directory had written the same five lines, already drifting on how
fieldnames are chosen and on what happens for an empty row list. Stdlib only on purpose:
`target_sites.py` pulls in gemmi and biopython, and a CSV pooler should not inherit a
structural-biology dependency to borrow a `DictWriter` call -- which also lets the generators
that run in the PyRosetta environment, where gemmi is absent, share this.

`read_designs` is here for the same reason. The name-to-sequence contract between
`collect_designs.py`, the histidine generators and `ph_switch.py` had been written twice and had
already drifted on how an empty filter is spelled; a design table read two ways is a funnel
whose denominators do not have to agree.

Named `csvio` rather than `tables` because `tables` is PyTables' import name on PyPI, and a
sibling-directory module that shadows an installed package is a trap worth not setting.
"""

from __future__ import annotations

import csv
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Any


def read_designs(path: Path, only: Collection[str] | None = None) -> dict[str, str]:
    """Name to sequence from a design table, restricted to `only` when it is non-empty."""
    with path.open(newline="") as handle:
        return {
            row["name"]: row["sequence"].strip().upper()
            for row in csv.DictReader(handle)
            if not only or row["name"] in only
        }


def write_rows(path: Path, rows: Sequence[Mapping[str, Any]],
               columns: Sequence[str] | None = None) -> None:
    """Write `rows` to a CSV, keyed on `columns` if given, else on the first row's keys.

    An empty `rows` writes a header-only file when `columns` says what the columns are, and
    an empty file otherwise -- rather than raising `IndexError` on `rows[0]`, which is how
    a stage with nothing to report used to fail.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(columns) if columns else list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
