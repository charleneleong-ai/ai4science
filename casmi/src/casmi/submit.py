"""Build a local CASMI submission CSV: 25 ranked SMILES per test molecule.

An earlier version also offered a "lookup" source — joining test spectra to bit-identical
train spectra and reading off their structures — on the belief that this recovered the
answers. It did not. Those structures are mostly not the competition's answers: a submission
whose top-1 matched them for 215 of 400 molecules scored 0.115, which bounds them correct
for at most ~46. The lookup is removed rather than left as a trap.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

from .harness import Query, Structure, retrieve

TOP_K = 25
TEST_COLUMNS = ["molecule_id", "adduct", "precursor_mz", "ms2_mzs", "ms2_normalized_intensities"]


@dataclass
class TestMolecule:
    molecule_id: str
    adduct: str
    precursor_mz: float
    spectra: list[tuple[list[float], list[float]]]

    def as_query(self) -> Query:
        q = Query(self.molecule_id, "", self.adduct, self.precursor_mz)
        q.spectra.extend(self.spectra)
        return q


def read_test(test: Path) -> list[TestMolecule]:
    """Group test spectra by molecule — the metric scores per molecule, not per spectrum."""
    table = pq.read_table(test, columns=TEST_COLUMNS)
    molecules: dict[str, TestMolecule] = {}
    for i in range(table.num_rows):
        mol_id = table.column("molecule_id")[i].as_py()
        adduct = table.column("adduct")[i].as_py()
        mz = table.column("precursor_mz")[i].as_py()
        mzs = table.column("ms2_mzs")[i].as_py()
        intens = table.column("ms2_normalized_intensities")[i].as_py()
        entry = molecules.setdefault(mol_id, TestMolecule(mol_id, adduct, mz, []))
        entry.spectra.append((mzs, intens))
    return list(molecules.values())


def pad(candidates: Sequence[str], k: int = TOP_K) -> list[str]:
    """Exactly k entries. MRR@k never penalises a wrong guess at a lower rank, so a short
    list only throws away free chances; repeat the last rather than leave slots empty."""
    out = [c for c in candidates[:k] if c]
    if not out:
        return ["C"] * k  # methane: a syntactically valid placeholder, never correct
    while len(out) < k:
        out.append(out[-1])
    return out


def rank_fallback(molecule: TestMolecule, pool: Sequence[Structure], masses: Sequence[float],
                  ranker, ppm: float) -> list[str]:
    query = molecule.as_query()
    mass = query.mass
    if mass is None:
        return []
    candidates = retrieve(pool, masses, mass, ppm)
    return [s.smiles for s in ranker(query, candidates)[:TOP_K]]


def write_submission(rows: dict[str, list[str]], out: Path) -> None:
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["molecule_id", "smiles"])
        for mol_id, smiles in rows.items():
            writer.writerow([mol_id, ";".join(smiles)])
