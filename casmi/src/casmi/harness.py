"""Molecule-disjoint evaluation for CASMI 2026.

The public leaderboard cannot be used for model selection: every public test spectrum is
duplicated in train (median peak-cosine 1.000 across all 400 molecules), so a score there
measures an exact join, not structure elucidation. See the 2026-10-06 learning-log entry.

So evaluation happens here instead, and it separates the two failure modes that a single
MRR number hides:

  recall  — is the true structure in the candidate pool at all? (a retrieval problem)
  ranking — given it is, does the scorer put it near the top? (the actual contest)

Splitting is on `inchikey14`, which is the connectivity block the competition metric itself
compares, so the split matches the equivalence class being scored.
"""

from __future__ import annotations

import bisect
import csv
import hashlib
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pyarrow.parquet as pq

from .adducts import neutral_mass, ppm_window
from .formula import formula_mass

QUERY_COLUMNS = ["inchikey14", "normalized_smiles", "adduct", "precursor_mz",
                 "ms2_mzs", "ms2_normalized_intensities", "instrument_type"]


@dataclass(frozen=True, slots=True)
class Structure:
    """A pool entry: a candidate the ranker may propose.

    Slotted because a realistic pool is ~740k of these and the per-object dict dominates.
    """

    inchikey14: str
    smiles: str
    mass: float


@dataclass
class Query:
    """One held-out molecule, with every spectrum recorded for it.

    The metric scores per molecule, not per spectrum, so aggregation across a molecule's
    spectra is part of the task rather than an optimisation.
    """

    inchikey14: str
    smiles: str
    adduct: str
    precursor_mz: float
    spectra: list[tuple[list[float], list[float]]] = field(default_factory=list)

    @property
    def mass(self) -> float | None:
        return neutral_mass(self.precursor_mz, self.adduct)


def holdout_fraction(inchikey14: str) -> float:
    """A stable per-molecule hash in [0, 1) — the split must not move between runs."""
    digest = hashlib.sha256(inchikey14.encode()).digest()
    return int.from_bytes(digest[:8], "big") / 2**64


def is_held_out(inchikey14: str, frac: float) -> bool:
    return holdout_fraction(inchikey14) < frac


def iter_batches(path: Path, columns: Sequence[str], batch_size: int = 100_000) -> Iterator:
    for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=list(columns)):
        yield batch


def build_pool(path: Path) -> list[Structure]:
    """Every unique structure in the file, as the candidate database.

    Masses come from `molecular_formula`, never from an observed precursor — see
    `casmi.formula` for why that distinction decides whether the harness measures anything.

    Deliberately includes held-out molecules: the realistic setting is that the structure
    exists in some database while no *spectrum* for it is available. A pool that excluded
    them would make recall 0 by construction and measure nothing.
    """
    seen: dict[str, Structure] = {}
    cols = ["inchikey14", "normalized_smiles", "molecular_formula"]
    for batch in iter_batches(path, cols):
        keys = batch.column("inchikey14").to_pylist()
        smiles = batch.column("normalized_smiles").to_pylist()
        formulas = batch.column("molecular_formula").to_pylist()
        for key, smi, formula in zip(keys, smiles, formulas, strict=True):
            if key is None or key in seen:
                continue
            mass = formula_mass(formula) if formula else None
            if mass is not None:
                seen[key] = Structure(key, smi, mass)
    return sorted(seen.values(), key=lambda s: s.mass)


def build_coconut_pool(csv_path: Path) -> list[Structure]:
    """COCONUT (CC0, ~740k natural products) as the candidate database.

    The train-only pool is not a fair test: it holds the answers because they are in the same
    file, and at 276k it is smaller than any database a real submission would search. This is
    the realism check — masses still come from `molecular_formula`, so the two pools are
    measured on one mass convention.
    """
    seen: dict[str, Structure] = {}
    with csv_path.open(newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            key = (row.get("standard_inchi_key") or "")[:14]
            smiles = row.get("canonical_smiles") or ""
            if len(key) != 14 or not smiles or key in seen:
                continue
            mass = formula_mass(row.get("molecular_formula") or "")
            if mass is not None:
                seen[key] = Structure(key, smiles, mass)
    return sorted(seen.values(), key=lambda s: s.mass)


def merge_pools(*pools: Sequence[Structure]) -> list[Structure]:
    """Union by inchikey14, mass-sorted. Earlier pools win on collision."""
    seen: dict[str, Structure] = {}
    for pool in pools:
        for s in pool:
            seen.setdefault(s.inchikey14, s)
    return sorted(seen.values(), key=lambda s: s.mass)


def collect_queries(path: Path, frac: float, limit: int = 0) -> list[Query]:
    """The held-out molecules and their spectra."""
    queries: dict[str, Query] = {}
    for batch in iter_batches(path, QUERY_COLUMNS):
        keys = batch.column("inchikey14").to_pylist()
        for i, key in enumerate(keys):
            if key is None or not is_held_out(key, frac):
                continue
            if key not in queries:
                if limit and len(queries) >= limit:
                    continue
                queries[key] = Query(
                    inchikey14=key,
                    smiles=batch.column("normalized_smiles")[i].as_py(),
                    adduct=batch.column("adduct")[i].as_py(),
                    precursor_mz=batch.column("precursor_mz")[i].as_py(),
                )
            queries[key].spectra.append(
                (batch.column("ms2_mzs")[i].as_py(), batch.column("ms2_normalized_intensities")[i].as_py())
            )
    return list(queries.values())


def retrieve(pool: Sequence[Structure], masses: Sequence[float], mass: float,
             ppm: float) -> Sequence[Structure]:
    """Pool entries inside the neutral-mass window.

    `pool` must be mass-sorted and `masses` its mass column — passed in rather than derived,
    so a per-query call doesn't rebuild a 275k-element list each time.
    """
    lo, hi = ppm_window(mass, ppm)
    return pool[bisect.bisect_left(masses, lo):bisect.bisect_right(masses, hi)]


Ranker = Callable[[Query, Sequence[Structure]], list[Structure]]


def rank_by_mass_error(query: Query, candidates: Sequence[Structure]) -> list[Structure]:
    """The floor: order by |mass error| alone, using no spectral information.

    Any real ranker must beat this, and it is not a trivial bar — the window is narrow, so
    on small candidate sets mass error alone can look deceptively strong.
    """
    mass = query.mass
    if mass is None:
        return list(candidates)
    return sorted(candidates, key=lambda s: abs(s.mass - mass))


def reciprocal_rank(ranked: Sequence[Structure], truth: str, k: int = 25) -> float:
    for i, s in enumerate(ranked[:k], start=1):
        if s.inchikey14 == truth:
            return 1.0 / i
    return 0.0


@dataclass
class Report:
    n_queries: int
    n_unmapped_adduct: int
    recall: float
    mrr: float
    median_candidates: float
    median_isomers: float

    def render(self) -> str:
        return (
            f"queries              {self.n_queries}\n"
            f"unmapped adduct      {self.n_unmapped_adduct}\n"
            f"candidate recall     {self.recall:.4f}   <- ceiling: ranking cannot exceed this\n"
            f"MRR@25               {self.mrr:.4f}\n"
            f"median candidates    {self.median_candidates:.0f}\n"
            f"median isomers       {self.median_isomers:.0f}   <- same formula as the truth, so\n"
            f"                             mass accuracy cannot separate them"
        )


def evaluate(queries: Sequence[Query], pool: Sequence[Structure], ppm: float = 5.0,
             ranker: Ranker = rank_by_mass_error, k: int = 25) -> Report:
    """Recall and MRR@k over held-out molecules, reported separately on purpose."""
    masses = [s.mass for s in pool]
    hits = 0
    rr_total = 0.0
    unmapped = 0
    counts: list[int] = []
    isomers: list[int] = []
    for query in queries:
        mass = query.mass
        if mass is None:
            unmapped += 1
            continue
        candidates = retrieve(pool, masses, mass, ppm)
        counts.append(len(candidates))
        truth = next((s for s in candidates if s.inchikey14 == query.inchikey14), None)
        if truth is None:
            continue
        hits += 1
        # candidates at the truth's exact formula mass: mass accuracy is blind between these,
        # so this is the ranking problem a spectral model actually has to solve
        isomers.append(sum(abs(s.mass - truth.mass) < 1e-6 for s in candidates))
        rr_total += reciprocal_rank(ranker(query, candidates), query.inchikey14, k)
    scored = len(counts) or 1
    counts.sort()
    isomers.sort()
    return Report(
        n_queries=len(queries),
        n_unmapped_adduct=unmapped,
        recall=hits / scored,
        mrr=rr_total / scored,
        median_candidates=counts[len(counts) // 2] if counts else 0.0,
        median_isomers=isomers[len(isomers) // 2] if isomers else 0.0,
    )
