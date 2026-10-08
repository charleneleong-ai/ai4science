"""CASMI 2026 Kaggle notebook — paste as a Script-type notebook and Submit.

Code competition, so: notebook only, **internet disabled**, writes
`/kaggle/working/submission.csv`. That constraint shapes the whole approach — no COCONUT or
PubChem fetch at inference time, so the candidate pool must be something already mounted.
Here it is `train.parquet`'s own 275k structures — and that is the flaw. **This notebook
scored 0.115 (rank 2397/2686).** The competition's answers are largely not train structures,
so a train-only pool cannot contain them. A real submission needs a PubChem-scale pool
attached as a Kaggle Dataset; the 0.421 public notebook does exactly that.

Self-contained on purpose. `pip install` needs internet, so the `casmi` package logic is
vendored here rather than imported; `tests/test_notebook_parity.py` in the repo pins this
copy against the package so the two cannot drift silently.

Approach, and why:
  pool      every unique structure in train, mass from `molecular_formula` — NEVER from an
            observed precursor m/z, which would make the true candidate float-identical to
            the query and score a measurement identity rather than chemistry.
  retrieve  5 ppm neutral-mass window. ~46 candidates, of which ~29 share the truth's exact
            formula, so mass accuracy is blind among them and ranking has to do the work.
  rank      fragment explainability: cut <= 2 bonds (ring bonds included — 68% of bonds in
            this chemistry are ring bonds), and score the share of peak intensity a
            structure's fragments can explain.

A "local estimate of 0.664" was quoted here before submission. It was wrong: it scored the
ranker against labels drawn from the same train pool the ranker searched, a closed loop with
no contact with ground truth. The same scorer rated the 0.421 public notebook at 0.030. The
mechanics below are sound and reproduce exactly on Kaggle; the pool is what fails.
"""

from __future__ import annotations

import bisect
import csv
import re
import subprocess
import sys
import time
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import pyarrow.parquet as pq


def ensure_rdkit() -> None:
    """Install rdkit from an attached wheel.

    rdkit is NOT in Kaggle's image — the first run died on `ModuleNotFoundError: rdkit` — and
    internet is off, so pip cannot reach PyPI. The wheel is attached as a dataset instead and
    installed with --no-index, which needs no network. cp313 because the worker runs 3.13.
    """
    try:
        import rdkit  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    wheels = sorted(Path("/kaggle/input").glob("**/rdkit-*.whl"))
    if not wheels:
        raise SystemExit(
            "rdkit is absent and no rdkit wheel is attached. Add the wheel dataset to this "
            "notebook's inputs — with internet off, pip cannot fetch it."
        )
    print(f"installing {wheels[0].name}", flush=True)
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-index", "--quiet", str(wheels[0])],
        check=True,
    )


ensure_rdkit()

from rdkit import Chem, RDLogger  # noqa: E402  (import follows the wheel install above)

RDLogger.DisableLog("rdApp.*")

PPM = 5.0
TOP_K = 25
MAX_BREAKS = 2
MAX_BONDS = 60
TOL_DA = 0.01
H_MASS = 1.00782503
H_SHIFTS = (-1, 0, 1)
PROTON = 1.007276

MONOISOTOPIC = {
    "H": 1.00782503, "D": 2.01410178, "C": 12.0, "N": 14.00307400, "O": 15.99491462,
    "F": 18.99840320, "Na": 22.98976928, "Si": 27.97692653, "P": 30.97376151,
    "S": 31.97207069, "Cl": 34.96885271, "K": 38.96370649, "Ca": 39.96259086,
    "Fe": 55.93494200, "Co": 58.93319400, "Ni": 57.93534241, "Cu": 62.92959772,
    "Zn": 63.92914201, "As": 74.92159460, "Se": 79.91651960, "Br": 78.91833760,
    "I": 126.90446800, "B": 11.00930550, "Mg": 23.98504170, "Mn": 54.93804391,
    "Hg": 201.97064340, "Pt": 194.96479170, "Sn": 119.90220160,
}
TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")

# adduct -> (multiplier on M, m/z offset); neutral M = (mz - offset) / multiplier
ADDUCTS = {
    "[M+H]+": (1, PROTON), "[M-H]-": (1, -PROTON), "[M+Na]+": (1, 22.989218),
    "[M+NH4]+": (1, 18.033823), "[M+K]+": (1, 38.963158),
    "[M+CH2O2-H]-": (1, 46.005480 - PROTON), "[M+Cl]-": (1, 34.969402),
    "[2M+H]+": (2, PROTON), "[2M+Na]+": (2, 22.989218),
}


def formula_mass(formula: str) -> float | None:
    """Monoisotopic mass, or None rather than a partial parse — a dropped element understates
    the mass and moves the candidate into the wrong retrieval window."""
    if not formula:
        return None
    cleaned = formula.strip().replace(" ", "").rstrip("+-")
    if not cleaned:
        return None
    total, pos = 0.0, 0
    for m in TOKEN.finditer(cleaned):
        if m.start() != pos:
            return None
        pos = m.end()
        element, digits = m.group(1), m.group(2)
        if element not in MONOISOTOPIC:
            return None
        total += MONOISOTOPIC[element] * (int(digits) if digits else 1)
    return total if pos == len(cleaned) and total else None


def neutral_mass(mz: float, adduct: str) -> float | None:
    entry = ADDUCTS.get(adduct)
    if entry is None:
        return None
    multiplier, offset = entry
    return (mz - offset) / multiplier


def component_masses(masses: list[float], bonds: list[tuple[int, int]], cut: set[int]) -> list[float]:
    """Piece masses after removing these bonds, by union-find over the remaining graph."""
    parent = list(range(len(masses)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, (a, b) in enumerate(bonds):
        if i in cut:
            continue
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    sums: dict[int, float] = defaultdict(float)
    for idx, mass in enumerate(masses):
        sums[find(idx)] += mass
    # round so differently-ordered sums don't persist as distinct floats and bloat the set
    return [round(v, 6) for v in sums.values()]


@lru_cache(maxsize=400_000)
def fragment_masses(smiles: str) -> tuple[float, ...]:
    """Sorted fragment masses reachable by cutting <= MAX_BREAKS bonds, plus the intact mass.

    Ring bonds are breakable and a ring needs two cuts to open, which is why MAX_BREAKS is 2:
    excluding rings left only ~19% of peak intensity explainable and the ranker dead.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ()
    atom_masses = [
        MONOISOTOPIC.get(a.GetSymbol(), 0.0) + a.GetTotalNumHs() * H_MASS for a in mol.GetAtoms()
    ]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    out = {round(sum(atom_masses), 6)}
    if len(bonds) > MAX_BONDS:
        return tuple(sorted(out))
    for n in range(1, MAX_BREAKS + 1):
        for combo in combinations(range(len(bonds)), n):
            out.update(component_masses(atom_masses, bonds, set(combo)))
    return tuple(sorted(out))


def any_within(sorted_masses: tuple[float, ...], target: float) -> bool:
    for shift in H_SHIFTS:
        want = target + shift * H_MASS
        i = bisect.bisect_left(sorted_masses, want - TOL_DA)
        if i < len(sorted_masses) and sorted_masses[i] <= want + TOL_DA:
            return True
    return False


def explained_fraction(mzs, intens, frags: tuple[float, ...], positive: bool) -> float:
    """Share of peak INTENSITY explained — one base peak says more than several noise peaks."""
    if not mzs or not frags:
        return 0.0
    charge = PROTON if positive else -PROTON
    total = sum(intens) or 1.0
    return sum(i for mz, i in zip(mzs, intens) if any_within(frags, mz - charge)) / total


def find_input() -> Path:
    """Kaggle mounts competition data under /kaggle/input/<slug>/."""
    for base in (Path("/kaggle/input"), Path.cwd()):
        hits = sorted(base.glob("**/train.parquet"))
        if hits:
            return hits[0].parent
    raise SystemExit("train.parquet not found under /kaggle/input")


def output_dir() -> Path:
    """`/kaggle/working` on Kaggle, the cwd anywhere else.

    Hardcoding the Kaggle path made this script impossible to run locally, which is how the
    first dry run died — after completing all the real work.
    """
    kaggle = Path("/kaggle/working")
    if kaggle.is_dir():
        return kaggle
    return Path.cwd()


def build_pool(train: Path) -> tuple[list[tuple[str, tuple[str, float]]], list[float]]:
    """Every unique structure in train, mass-sorted, with its mass column.

    Masses come from `molecular_formula`, never an observed precursor m/z — the latter makes
    the true candidate float-identical to the query and scores measurement identity.
    """
    pool: dict[str, tuple[str, float]] = {}
    for batch in pq.ParquetFile(train).iter_batches(
            batch_size=200_000, columns=["inchikey14", "normalized_smiles", "molecular_formula"]):
        keys = batch.column("inchikey14").to_pylist()
        smis = batch.column("normalized_smiles").to_pylist()
        forms = batch.column("molecular_formula").to_pylist()
        for key, smi, form in zip(keys, smis, forms):
            if not key or key in pool or not smi:
                continue
            mass = formula_mass(form or "")
            if mass is not None:
                pool[key] = (smi, mass)
    entries = sorted(pool.items(), key=lambda kv: kv[1][1])
    return entries, [v[1] for _, v in entries]


def read_test(test: Path) -> dict[str, dict]:
    """Test spectra grouped by molecule — the metric scores per molecule, not per spectrum."""
    t = pq.read_table(test, columns=["molecule_id", "adduct", "precursor_mz",
                                     "ms2_mzs", "ms2_normalized_intensities"])
    molecules: dict[str, dict] = {}
    for i in range(t.num_rows):
        mol_id = t.column("molecule_id")[i].as_py()
        entry = molecules.setdefault(mol_id, {"adduct": t.column("adduct")[i].as_py(),
                                              "mz": t.column("precursor_mz")[i].as_py(),
                                              "spectra": []})
        entry["spectra"].append((t.column("ms2_mzs")[i].as_py(),
                                 t.column("ms2_normalized_intensities")[i].as_py()))
    return molecules


def rank_one(mol: dict, entries, pool_masses: list[float]) -> list[str]:
    """Retrieve by neutral-mass window, then order by explained intensity."""
    mass = neutral_mass(mol["mz"], mol["adduct"]) if mol["mz"] else None
    if mass is None:
        return []
    tol = mass * PPM * 1e-6
    lo = bisect.bisect_left(pool_masses, mass - tol)
    hi = bisect.bisect_right(pool_masses, mass + tol)
    positive = mol["adduct"].endswith("+")
    scored = []
    for _key, (smi, cand_mass) in entries[lo:hi]:
        frags = fragment_masses(smi)
        best = max((explained_fraction(mzs, ins, frags, positive)
                    for mzs, ins in mol["spectra"]), default=0.0)
        scored.append((-best, abs(cand_mass - mass), smi))
    # Sort on the score keys ONLY. A bare .sort() falls through to the SMILES string, breaking
    # ties alphabetically — and with ~29 isomers sharing an identical formula mass, ties are
    # the common case, so that silently diverged from the package (202/400 rows) and
    # made the notebook submit something other than what was measured. A stable sort on (-score, mass error) preserves the mass
    # order the slice already has, matching casmi.fragments.rank_by_fragments exactly.
    scored.sort(key=lambda row: (row[0], row[1]))
    return [smi for _, _, smi in scored[:TOP_K]]


def pad(picks: list[str]) -> list[str]:
    """Exactly TOP_K non-empty entries — a null or short row is rejected outright, and a
    wrong guess at a lower rank costs nothing in MRR@k."""
    out = [p for p in picks[:TOP_K] if p] or ["C"]
    while len(out) < TOP_K:
        out.append(out[-1])
    return out


def write_submission(rows: dict[str, list[str]], expected: int) -> Path:
    """Validate against the stated rejection rules, then write."""
    assert rows, "no rows produced"
    assert len(rows) == expected, "molecule_id count changed"
    for mol_id, picks in rows.items():
        assert mol_id and len(picks) == TOP_K, f"bad row for {mol_id}"
        assert all(p for p in picks), f"empty candidate for {mol_id}"
    out = output_dir() / "submission.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["molecule_id", "smiles"])
        for mol_id, picks in rows.items():
            w.writerow([mol_id, ";".join(picks)])
    return out


def main() -> None:
    t0 = time.time()
    data = find_input()
    print(f"input: {data}", flush=True)

    entries, pool_masses = build_pool(data / "train.parquet")
    print(f"pool: {len(entries):,} structures  ({time.time()-t0:.0f}s)", flush=True)

    molecules = read_test(data / "test.parquet")
    print(f"test molecules: {len(molecules)}", flush=True)

    rows: dict[str, list[str]] = {}
    for n, (mol_id, mol) in enumerate(molecules.items(), 1):
        rows[mol_id] = pad(rank_one(mol, entries, pool_masses))
        if n % 50 == 0:
            print(f"  ranked {n}/{len(molecules)}  ({time.time()-t0:.0f}s)", flush=True)

    out = write_submission(rows, len(molecules))
    print(f"wrote {out}: {len(rows)} rows x {TOP_K}  ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
