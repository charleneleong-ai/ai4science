"""CASMI 2026 Kaggle notebook. Push with `kaggle kernels push -p casmi/notebooks`; the inputs it
needs are listed in `kernel-metadata.json` beside it.

Code competition, so: notebook only, **internet disabled**, writes
`/kaggle/working/submission.csv`. Every pool must therefore be mounted as a dataset.

Scores so far, each changing one thing:
  0.115  train's 275k structures only — the answers are largely not train structures
  0.049  + the 106M-structure PubChem tier, 4000 nearest-mass per molecule — arbitrary among
         same-formula isomers, and fragments alone cannot find the answer among thousands
  0.135  PubChem shortlist chosen by a popularity prior, ranked by fragments + prior
  0.250  the same, with the public FPNet's fingerprint score in place of fragments
  0.258  the same, averaging two FPNet checkpoints (full1 + FPNet A)
  0.254    with popularity weight 0.15 instead of 0.25 — kept at 0.25
  0.261    with a 1000-connectivity shortlist instead of 300 — adopted
  0.262  the same, never proposing the train structure behind a copied spectrum
  0.286  the same, with ICEBERG re-ordering same-formula groups in the top 60 (CPU: 269 / 400
         molecules scored inside the 6 h budget)
  0.292  the same on a T4: all 371 coverable molecules scored by ICEBERG
  0.320  with GLACIER fused beside ICEBERG

Self-contained on purpose. `pip install` needs internet, so the `casmi` package logic is
vendored here rather than imported; `tests/test_notebook_parity.py` in the repo pins this
copy against the package so the two cannot drift silently.

Approach:
  pool      train (mass from `molecular_formula`, never an observed precursor m/z) plus the
            SHORTLIST most-documented PubChem connectivities inside 5 ppm.
  rank      z-scored spectral score plus POP_WEIGHT x log1p(substances) + log1p(PubMed
            articles). The spectral score is FPNet's fingerprint match when its datasets are
            attached, else fragment explainability (<= 2 bond cuts, ring bonds included).
            Compounds that reach a mass spectrometer are usually ones that were bought,
            isolated and written about.

The local scorer once quoted 0.664 for a submission that scored 0.115; it was never calibrated
and nothing here is scored locally. Only the leaderboard measures this notebook.
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
from typing import NamedTuple

import numpy as np
import pyarrow.parquet as pq


def ensure_rdkit() -> None:
    """Install rdkit from an attached wheel.

    rdkit is NOT in Kaggle's image — the first run died on `ModuleNotFoundError: rdkit` — and
    internet is off, so pip cannot reach PyPI. The wheel is attached as a dataset instead and
    installed with --no-index, which needs no network. The wheel must match the interpreter: the
    pinned image (kernel-metadata.json) runs the 3.12 that ICEBERG's own wheels need.
    """
    try:
        import rdkit  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
    wheels = sorted(Path("/kaggle/input").glob(f"**/rdkit-*-{tag}-{tag}-*.whl"))
    if not wheels:
        raise SystemExit(
            f"rdkit is absent and no {tag} rdkit wheel is attached. Add the wheel dataset to "
            "this notebook's inputs — with internet off, pip cannot fetch it."
        )
    print(f"installing {wheels[0].name}", flush=True)
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-index", "--quiet", str(wheels[0])],
        check=True,
    )


ensure_rdkit()

from rdkit import Chem, RDLogger  # noqa: E402  (import follows the wheel install above)
from rdkit.Chem.Descriptors import ExactMolWt  # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PPM = 5.0
TOP_K = 25
SHORTLIST = 1000  # distinct connectivities per molecule, chosen by popularity; 300 scored 0.258
RERANK_N = 60  # ranked list depth handed to ICEBERG; the dataset README's recommendation
ICE_LAMBDA = 0.5  # weights on each forward model's z-score within a same-formula group, as recommended
GL_LAMBDA = 0.5
ICE_BUDGET_S = 6 * 3600  # caps: ICEBERG took 20 min on a T4 and covered 269 / 371 in 6 h on CPU;
GL_BUDGET_S = 3 * 3600  # the runners stop cleanly at them. GLACIER runs only on a GPU.
POP_WEIGHT = 0.25  # on a z-scored spectral score, as in the published prior's recipe
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
    """Test spectra grouped by molecule — the metric scores per molecule, not per spectrum.

    `records` carries what FPNet conditions on (precursor, adduct, polarity, collision energy),
    in the shape the public pipeline builds them.
    """
    cols = ["molecule_id", "adduct", "precursor_mz", "ms2_mzs", "ms2_normalized_intensities",
            "ionization_mode", "collision_energy_ev"]
    molecules: dict[str, dict] = {}
    for r in pq.read_table(test, columns=cols).to_pylist():
        entry = molecules.setdefault(r["molecule_id"], {"adduct": r["adduct"], "mz": r["precursor_mz"],
                                                        "records": []})
        ces = r["collision_energy_ev"] or []
        entry["records"].append(dict(
            mz=r["ms2_mzs"], it=r["ms2_normalized_intensities"], prec=float(r["precursor_mz"]),
            adduct=r["adduct"], mode=1 if r["ionization_mode"] == "positive" else -1,
            ce=float(np.mean(ces)) if ces else None, ce_n=len(ces) or 1))
    return molecules


class Candidate(NamedTuple):
    smiles: str
    mass: float
    pop: float


class PubChemTier:
    """The community `casmi26-pubchem-tier` (~106M structures, mass-sorted) with its popularity prior.

    Memory-mapped so 7.2 GB never has to fit in RAM: `pc_mass` (sorted), `pc_off` (N+1 byte
    offsets), `pc_smiles` (one flat byte buffer). `casmi26-pubchem-popularity-prior` adds row-aligned
    `pc_lsid` / `pc_lpmid` (log1p of PubChem substance records / PubMed articles) and `pc_ik14`.

    A 5 ppm window holds more than 4000 structures for 297/400 molecules. Cutting it by mass
    error is arbitrary among same-formula isomers and scored 0.049; the prior picks instead.
    """

    def __init__(self, root: Path, prior: Path) -> None:
        self.mass = np.load(root / "pc_mass.npy", mmap_mode="r")
        self.off = np.load(root / "pc_off.npy", mmap_mode="r")
        self.smiles = np.load(root / "pc_smiles.npy", mmap_mode="r")
        self.lsid = np.load(prior / "pc_lsid.npy", mmap_mode="r")
        self.lpmid = np.load(prior / "pc_lpmid.npy", mmap_mode="r")
        self.ik14 = np.load(prior / "pc_ik14.npy", mmap_mode="r")
        assert len(self.off) == len(self.mass) + 1, "expected N+1 offsets for N structures"
        assert len(self.lsid) == len(self.lpmid) == len(self.ik14) == len(self.mass), \
            "prior arrays are not row-aligned with the tier"

    def smiles_at(self, i: int) -> str:
        return bytes(self.smiles[int(self.off[i]):int(self.off[i + 1])]).decode()

    def window(self, lo: float, hi: float,
               train_keys: list[str]) -> tuple[dict[str, Candidate], dict[str, float]]:
        """The SHORTLIST most popular connectivities in [lo, hi], and the popularity of each
        train key found anywhere in the window — not only in the shortlist, or every train
        structure outside it would be ranked as if PubChem had never heard of it."""
        a = int(np.searchsorted(self.mass, lo, "left"))
        b = int(np.searchsorted(self.mass, hi, "right"))
        keys = np.asarray(self.ik14[a:b])
        masses = self.mass[a:b]
        pop = np.add(self.lsid[a:b], self.lpmid[a:b], dtype=np.float32)
        shortlist = {
            keys[j].decode(): Candidate(self.smiles_at(a + j), float(masses[j]), float(pop[j]))
            for j in top_distinct(pop, keys, SHORTLIST)
        }
        train_pop: dict[str, float] = {}
        want = np.array([k.encode() for k in train_keys], dtype="S14")
        for j in np.flatnonzero(np.isin(keys, want)):
            key = keys[j].decode()
            train_pop[key] = max(train_pop.get(key, 0.0), float(pop[j]))
        return shortlist, train_pop

    def mass_convention_error(self, samples: int = 60) -> float:
        """Median |pc_mass - RDKit exact mass| over a spread sample, in Da.

        If pc_mass were average molecular weight, every ppm window would land in the wrong
        place and retrieval would silently miss, so this is checked rather than assumed.
        """
        errors = []
        for i in np.linspace(0, len(self.mass) - 1, samples, dtype=np.int64):
            mol = Chem.MolFromSmiles(self.smiles_at(int(i)))
            if mol is not None:
                errors.append(abs(float(self.mass[i]) - ExactMolWt(mol)))
        return float(np.median(errors)) if errors else float("inf")


def top_distinct(pop: np.ndarray, keys: np.ndarray, k: int) -> list[int]:
    """Indices of the k most popular rows with distinct keys; the first (most popular) row wins.

    The metric scores InChIKey14, so stereoisomers and tautomers sharing one would only take
    each other's slots. Rows without a key have no prior and are left out.
    """
    picked: list[int] = []
    seen: set[bytes] = {b""}
    for j in np.argsort(-pop, kind="stable"):
        if keys[j] in seen:
            continue
        seen.add(keys[j])
        picked.append(int(j))
        if len(picked) == k:
            break
    return picked


def find_pubchem() -> PubChemTier | None:
    tier = sorted(Path("/kaggle/input").glob("**/pc_mass.npy"))
    prior = sorted(Path("/kaggle/input").glob("**/pc_lsid.npy"))
    if not tier:
        return None
    if not prior:
        raise SystemExit("PubChem tier attached without casmi26-pubchem-popularity-prior.")
    return PubChemTier(tier[0].parent, prior[0].parent)


class FPNetScorer:
    """Spectrum -> fingerprint logits z (the public FPNet); a candidate with fingerprint f scores f.z.

    Model code is the public `casmi` package in casmi26-v4b-models; the checkpoints are
    casmi26-fpnet-full1 and that package's FPNet A, whose logits ModelBank averages (the 0.421
    notebook's 'ens' bank); the 10,226 informative fingerprint bits
    come from casmi26-v2-pool. Single-spectrum and per-polarity merged logits are averaged, as in
    that package's Engine. Imported in __init__, not at the top: torch, numba and that package exist
    only where the datasets are mounted.
    """

    def __init__(self, code: Path, ckpts: list[Path], bits: Path) -> None:
        sys.path.insert(0, str(code))
        import torch
        from casmi import chem, fpnet
        from casmi.spectra import merge_spectra

        self.chem, self.fpnet, self.merge = chem, fpnet, merge_spectra
        self.bank = fpnet.ModelBank([str(c) for c in ckpts],
                                    device="cuda" if torch.cuda.is_available() else "cpu")
        self.bits = np.load(bits)
        assert len(self.bits) == self.bank.nbits, "fp_bits does not match the checkpoint's output"
        self._fp: dict[str, np.ndarray] = {}  # SMILES -> packed selected bits

    def item(self, mz, it, prec: float, adduct: str, ce: list[float], n: int, mode: float) -> dict:
        pm, pi = self.fpnet.prep_peaks(np.asarray(mz, np.float64), np.asarray(it, np.float64), prec)
        return dict(mz=pm, it=pi, prec=prec, adduct_ix=self.chem.adduct_index(adduct),
                    ce=float(np.mean(ce)) if ce else 0.0, ce_known=1.0 if ce else 0.0,
                    n_merged=min(n, 8), mode=mode)

    def logits(self, records: list[dict]) -> np.ndarray:
        single = [self.item(r["mz"], r["it"], r["prec"], r["adduct"],
                            [r["ce"]] if r["ce"] is not None else [], r["ce_n"], float(r["mode"]))
                  for r in records]
        merged = []
        for mode in (1, -1):
            grp = [r for r in records if r["mode"] == mode]
            if not grp:
                continue
            mz, it = self.merge([(r["mz"], r["it"]) for r in grp])
            adducts = [r["adduct"] for r in grp]
            merged.append(self.item(mz, it, float(np.median([r["prec"] for r in grp])),
                                    max(set(adducts), key=adducts.count),
                                    [r["ce"] for r in grp if r["ce"] is not None],
                                    sum(r["ce_n"] for r in grp), float(mode)))
        z = self.bank.logits(single + merged)
        return 0.5 * (z[:len(single)].mean(0) + z[len(single):].mean(0))

    def fingerprint(self, smi: str) -> np.ndarray:
        if smi not in self._fp:
            fp = self.chem.raw_fingerprint(smi)
            self._fp[smi] = np.packbits(fp[self.bits] if fp is not None else np.zeros(len(self.bits), np.uint8))
        return self._fp[smi]

    def fingerprints(self, smiles: list[str]) -> np.ndarray:
        packed = np.stack([self.fingerprint(s) for s in smiles])
        return np.unpackbits(packed, axis=1)[:, :len(self.bits)].astype(np.float32)

    def scores(self, records: list[dict], smiles: list[str]) -> np.ndarray:
        return self.fingerprints(smiles) @ self.logits(records)


def find_fpnet() -> FPNetScorer | None:
    """All the pieces or none: a partial attachment would silently fall back to fragments."""
    root = Path("/kaggle/input")
    found = [sorted(root.glob(pat)) for pat in
             ("**/code/casmi/fpnet.py", "**/fpnet_full1.pt", "**/casmi26-v4b-models/**/fpnet_0.pt",
              "**/casmi26-v2-pool/**/fp_bits.npy")]
    if not any(found):
        return None
    if not all(found):
        raise SystemExit("FPNet partly attached: need casmi26-v4b-models, casmi26-fpnet-full1 and "
                         "casmi26-v2-pool together.")
    return FPNetScorer(found[0][0].parent.parent, [found[1][0], found[2][0]], found[3][0])


def copied_parents(train: Path, molecules: dict[str, dict]) -> dict[str, set[str]]:
    """molecule_id -> InChIKey14s of the train rows whose spectrum is bit-identical to one of its own.

    1189 / 1213 test spectra are copies of train rows, but the train structure attached is the
    answer for at most ~46 of 400 molecules. FPNet full1 was trained on those very rows, so it
    returns the copied structure with confidence: the 0.261 submission's top-1 was the copied
    structure for 267 / 400, while the 0.421 public notebook's was for 12 / 400.
    """
    want: dict[tuple, set[str]] = {}
    for mol_id, mol in molecules.items():
        for r in mol["records"]:
            key = (r["adduct"], round(r["prec"] * 1e4), tuple(r["mz"]), tuple(r["it"]))
            want.setdefault(key, set()).add(mol_id)
    coarse = {k[:2] for k in want}
    out: dict[str, set[str]] = defaultdict(set)
    cols = ["inchikey14", "adduct", "precursor_mz", "ms2_mzs", "ms2_normalized_intensities"]
    for batch in pq.ParquetFile(train).iter_batches(batch_size=200_000, columns=cols):
        adducts = batch.column("adduct").to_pylist()
        precs = batch.column("precursor_mz").to_pylist()
        for i, (adduct, prec) in enumerate(zip(adducts, precs)):
            if prec is None or (adduct, round(prec * 1e4)) not in coarse:
                continue
            key = (adduct, round(prec * 1e4), tuple(batch.column("ms2_mzs")[i].as_py()),
                   tuple(batch.column("ms2_normalized_intensities")[i].as_py()))
            for mol_id in want.get(key, ()):
                out[mol_id].add(batch.column("inchikey14")[i].as_py())
    return dict(out)


def candidates_for(mass: float, entries, pool_masses: list[float], pubchem: PubChemTier | None,
                   exclude: set[str] = frozenset()) -> dict[str, Candidate]:
    """ik14 -> Candidate: the prior's PubChem shortlist plus the train window, minus `exclude`.

    Train structures keep their own SMILES and take their popularity from PubChem (0 if absent).
    """
    tol = mass * PPM * 1e-6
    lo = bisect.bisect_left(pool_masses, mass - tol)
    hi = bisect.bisect_right(pool_masses, mass + tol)
    train = entries[lo:hi]
    out: dict[str, Candidate] = {}
    train_pop: dict[str, float] = {}
    if pubchem is not None:
        out, train_pop = pubchem.window(mass - tol, mass + tol, [key for key, _ in train])
    for key, (smi, m) in train:
        out[key] = Candidate(smi, m, train_pop.get(key, 0.0))
    return {k: c for k, c in out.items() if k not in exclude}


def blend(spectral: np.ndarray, pop: np.ndarray) -> np.ndarray:
    """z-scored spectral score plus POP_WEIGHT * raw popularity, the published recipe.

    Popularity is not z-scored, so across a shortlist spanning ~5-20 it dominates and the spectral
    score mostly reorders near-ties. That is the hypothesis under test, not a tuned balance.
    """
    sd = spectral.std()
    z = (spectral - spectral.mean()) / sd if sd > 0 else np.zeros_like(spectral)
    return z + POP_WEIGHT * pop


def explained_scores(mol: dict, smiles: list[str]) -> np.ndarray:
    positive = mol["adduct"].endswith("+")
    out = np.zeros(len(smiles))
    for i, smi in enumerate(smiles):
        frags = fragment_masses(smi)
        out[i] = max((explained_fraction(r["mz"], r["it"], frags, positive)
                      for r in mol["records"]), default=0.0)
    return out


def rank_candidates(mol: dict, mass: float, candidates: dict[str, Candidate],
                    fpnet: FPNetScorer | None = None) -> list[tuple[str, float]]:
    """Order by blended score — FPNet's f.z when attached, else fragment explainability — with
    mass error breaking ties."""
    rows = list(candidates.values())
    smiles = [c.smiles for c in rows]
    spectral = fpnet.scores(mol["records"], smiles) if fpnet else explained_scores(mol, smiles)
    score = blend(spectral, np.array([c.pop for c in rows]))
    order = sorted(range(len(rows)), key=lambda i: (-score[i], abs(rows[i].mass - mass)))
    return [(rows[i].smiles, float(score[i])) for i in order[:RERANK_N]]


class ForwardRescorer:
    """Forward models re-order same-formula groups inside the top RERANK_N: ICEBERG (casmi26-iceberg)
    and, when attached, GLACIER (casmi26-glacier). Each predicts every candidate's spectrum and scores
    it against the query; the fusion is z(score) + ICE_LAMBDA z(ICEBERG) + GL_LAMBDA z(GLACIER)
    within each group, the datasets' recommended recipe.

    Both run in subprocesses with ICEBERG's wheels, including a cp312 RDKit, and return null scores
    on any failure without raising, so the logged counts are the only evidence they worked.
    Imported in __init__: the packages exist only on Kaggle.
    """

    def __init__(self, ice_pkg: Path, gl_pkg: Path | None) -> None:
        sys.path.append(str(ice_pkg))
        import fuse

        self.ice_pkg, self.gl_pkg, self.fuse, self.gl_fuse = ice_pkg, gl_pkg, fuse, None
        if gl_pkg:
            sys.path.append(str(gl_pkg))
            import gl_fuse

            self.gl_fuse = gl_fuse

    @property
    def names(self) -> str:
        return "ICEBERG + GLACIER" if self.gl_fuse else "ICEBERG"

    @staticmethod
    def formula(smi: str) -> str:  # an unparsable SMILES becomes its own group, so rerank never moves it
        mol = Chem.MolFromSmiles(smi)
        return CalcMolFormula(mol) if mol else smi

    def items(self, test: Path, ranked: dict[str, list[tuple[str, float]]],
              formulas: dict[str, list[str]]) -> list[dict]:
        """Runner input: only groups with a member in the top TOP_K can change the submission."""
        import pandas as pd

        cands = {}
        for m, r in ranked.items():
            smis = [s for s, _ in r]
            reachable = set(formulas[m][:TOP_K])
            group = set(self.fuse.ice_candidates(smis, formulas[m], RERANK_N))
            cands[m] = [s for s, f in zip(smis, formulas[m]) if s in group and f in reachable]
        return self.fuse.build_ice_input(pd.read_parquet(test), {m: c for m, c in cands.items() if c})

    def score(self, items: list[dict]) -> list[dict[str, dict[str, float | None]]]:
        """Scores per model that returned any; a model returning none is reported and left out of the fusion."""
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"forward models: {sum(bool(i['spectra']) for i in items)} molecules with covered spectra", flush=True)
        runs = [("ICEBERG", self.fuse.run_ice(str(self.ice_pkg), items, workdir="/tmp/ice_work", device=device,
                                              budget_s=ICE_BUDGET_S, site="/tmp/ice_site"))]
        if self.gl_fuse and device == "cuda":
            runs.append(("GLACIER", self.gl_fuse.run_gl(str(self.gl_pkg), items, workdir="/tmp/gl_work",
                                                        device=device, budget_s=GL_BUDGET_S, site="/tmp/ice_site",
                                                        wheels=str(self.ice_pkg / "wheels"))))
        kept = []
        for name, scores in runs:
            returned = sum(v is not None for per in scores.values() for v in per.values())
            print(f"{name}: {returned} candidate scores returned{'' if returned else ' — FAILED, left out'}",
                  flush=True)
            if returned:
                kept.append(scores)
        return kept

    def reorder(self, ranked: dict[str, list[tuple[str, float]]], formulas: dict[str, list[str]],
                runs: list[dict[str, dict[str, float | None]]]) -> dict[str, list[str]]:
        lams = [ICE_LAMBDA, GL_LAMBDA][:len(runs)]
        out = {}
        for m, r in ranked.items():
            smis = [s for s, _ in r]
            dicts = [run.get(m, {}) for run in runs] or [{}]
            order = (self.gl_fuse.rerank_multi(smis, None, [sc for _, sc in r], formulas[m], dicts, lams, top_n=RERANK_N)
                     if len(runs) > 1 else
                     self.fuse.rerank(smis, None, [sc for _, sc in r], formulas[m], dicts[0], lam=lams[0], top_n=RERANK_N))
            out[m] = [smis[i] for i in order]
        return out

    def rerank(self, test: Path, ranked: dict[str, list[tuple[str, float]]]) -> dict[str, list[str]]:
        formulas = {m: [self.formula(s) for s, _ in r] for m, r in ranked.items()}
        runs = self.score(self.items(test, ranked, formulas))
        if len(runs) > 1:  # ICEBERG alone, for comparison with the ICEBERG-only runs
            write_submission(self.reorder(ranked, formulas, runs[:1]), len(ranked), "iceberg_only.csv")
        out = self.reorder(ranked, formulas, runs)
        print(f"{self.names} changed the top-1 for "
              f"{sum(out[m][:1] != [s for s, _ in ranked[m][:1]] for m in ranked)} molecules", flush=True)
        return out


def find_forward() -> ForwardRescorer | None:
    ice = sorted(Path("/kaggle/input").glob("**/casmi26-iceberg/**/ice_runner.py"))
    gl = sorted(Path("/kaggle/input").glob("**/casmi26-glacier/**/gl_runner.py"))
    if gl and not ice:
        raise SystemExit("GLACIER attached without casmi26-iceberg, whose wheels and glue it needs.")
    return ForwardRescorer(ice[0].parent, gl[0].parent if gl else None) if ice else None


def pad(picks: list[str]) -> list[str]:
    """Exactly TOP_K non-empty entries — a null or short row is rejected outright, and a
    wrong guess at a lower rank costs nothing in MRR@k."""
    out = [p for p in picks[:TOP_K] if p] or ["C"]
    while len(out) < TOP_K:
        out.append(out[-1])
    return out


def write_submission(rows: dict[str, list[str]], expected: int, name: str = "submission.csv") -> Path:
    """Pad each row to TOP_K, validate against the stated rejection rules, then write."""
    rows = {m: pad(p) for m, p in rows.items()}
    assert rows, "no rows produced"
    assert len(rows) == expected, "molecule_id count changed"
    for mol_id, picks in rows.items():
        assert mol_id and len(picks) == TOP_K, f"bad row for {mol_id}"
        assert all(p for p in picks), f"empty candidate for {mol_id}"
    out = output_dir() / name
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["molecule_id", "smiles"])
        for mol_id, picks in rows.items():
            w.writerow([mol_id, ";".join(picks)])
    return out


def open_pubchem() -> PubChemTier | None:
    """Attach the PubChem tier if present, refusing to run on a mass convention that is wrong."""
    pubchem = find_pubchem()
    if pubchem is None:
        print("PubChem tier NOT attached — train-only pool, which scored 0.115", flush=True)
        return None
    err = pubchem.mass_convention_error()
    print(f"pubchem: {len(pubchem.mass):,} structures; median |pc_mass - exact| = {err:.6f} Da",
          flush=True)
    if err > 0.005:
        raise SystemExit(f"pc_mass is not monoisotopic (median error {err:.4f} Da); every ppm "
                         "window would miss. Stopping rather than ranking the wrong candidates.")
    return pubchem


def rank_all(molecules: dict[str, dict], entries, pool_masses: list[float], pubchem: PubChemTier | None,
             fpnet: FPNetScorer | None, copied: dict[str, set[str]]) -> dict[str, list[tuple[str, float]]]:
    t0 = time.time()
    ranked: dict[str, list[tuple[str, float]]] = {}
    sizes: list[int] = []
    for n, (mol_id, mol) in enumerate(molecules.items(), 1):
        mass = neutral_mass(mol["mz"], mol["adduct"]) if mol["mz"] else None
        ranked[mol_id] = []
        if mass is not None:
            cands = candidates_for(mass, entries, pool_masses, pubchem, copied.get(mol_id, set()))
            sizes.append(len(cands))
            ranked[mol_id] = rank_candidates(mol, mass, cands, fpnet)
        if n % 25 == 0:
            print(f"  ranked {n}/{len(molecules)}  median candidates {int(np.median(sizes))}  "
                  f"({time.time()-t0:.0f}s)", flush=True)
    print(f"candidates: median {int(np.median(sizes))}, max {max(sizes)}", flush=True)
    return ranked


def main() -> None:
    t0 = time.time()
    data = find_input()
    print(f"input: {data}", flush=True)

    entries, pool_masses = build_pool(data / "train.parquet")
    print(f"train pool: {len(entries):,} structures  ({time.time()-t0:.0f}s)", flush=True)
    pubchem = open_pubchem()
    fpnet = find_fpnet()
    print(f"spectral score: {'FPNet' if fpnet else 'fragment explainability'}", flush=True)
    forward = find_forward()
    print(f"forward model: {forward.names if forward else 'none'}", flush=True)

    molecules = read_test(data / "test.parquet")
    print(f"test molecules: {len(molecules)}", flush=True)
    copied = copied_parents(data / "train.parquet", molecules)
    print(f"copied train structures excluded for {len(copied)} molecules "
          f"({sum(map(len, copied.values()))} in all)  ({time.time()-t0:.0f}s)", flush=True)

    ranked = rank_all(molecules, entries, pool_masses, pubchem, fpnet, copied)
    final = {m: [s for s, _ in r] for m, r in ranked.items()}
    if forward:
        write_submission(final, len(molecules), "pre_forward.csv")
        final = forward.rerank(data / "test.parquet", ranked)
    out = write_submission(final, len(molecules))
    print(f"wrote {out}: {len(final)} rows x {TOP_K}  ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
