"""The table readers, writers and name contract this directory shares, in one place.

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
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any, TypeVar

T = TypeVar("T")
K = TypeVar("K")


def read_rows(path: Path) -> list[dict[str, str]]:
    """Every row of a CSV as a dict. Raises if the file is absent, which is the honest answer.

    `tnf_binder.tables.read_rows` returns an empty table for a missing file, on the grounds that
    a stage which has not run yet is the normal case. That is right for a pipeline that polls
    its own intermediate outputs and wrong here: every path this reads is a required CLI
    argument naming a completed measurement, and an empty table would silently become an empty
    funnel -- a selection made from no candidates, reported as a selection.
    """
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: str | None, default: float = 0.0) -> float:
    """A numeric cell, falling back to `default` when it is blank or malformed.

    BoltzGen, PROPKA and pKAI all leave cells empty rather than writing a zero, and a funnel
    that raises on the first empty cell reports nothing instead of reporting what it measured.
    """
    if not value:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def group_by(rows: Iterable[T], key: Callable[[T], K]) -> dict[K, list[T]]:
    """Rows bucketed by `key`, preserving the order within each bucket.

    `ph_switch.group_by` is the same three lines. It is not imported from there because that
    module pulls in gemmi, which the PyRosetta environment lacks; this is the home the rest of
    the directory can reach. TODO: have `ph_switch` import this one at the next full re-run of
    the pKa pass, rather than editing the module that produced the committed measurements now.
    """
    grouped: dict[K, list[T]] = {}
    for row in rows:
        grouped.setdefault(key(row), []).append(row)
    return grouped


def same_design(design: str, query: str) -> bool:
    """Whether `query` names `design`, when two stages spell one design differently.

    `collect_designs.py` prefixes each design with the spec directory that produced it while
    BoltzGen's own id already carries that prefix, so a design written to disk as
    `conserved_s_45.cif` is recorded in the pooled table as `conserved_s_conserved_s_45`. Both
    spellings have to resolve to one design, in either direction, and padding with underscores
    is what stops `conserved_s_4` claiming `conserved_s_45`.

    This is the third spelling of that one rule -- `ph_switch.index_structures` tests the name
    inside the stem, `denovo_his.structure_paths` tests the stem inside the name -- so it lives
    here. TODO: the root fix is to stop doubling the prefix at `collect_designs.py:155`
    (`if not row[identity].startswith(directory.name)`), after which `stem == design` and all
    three callers can be deleted. It is deliberately not done now: renaming every design
    invalidates the m8 novelty queries, the `design` keys in the agreement table and every
    variant name, which means re-running the measurements to delete a predicate.
    """
    stem = Path(query).stem
    return stem == design or f"_{stem}_" in f"_{design}_"


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
