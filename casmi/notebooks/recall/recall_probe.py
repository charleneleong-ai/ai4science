"""Kaggle kernel: how much a spectrum-driven second PubChem pass adds to the popularity shortlist.

The 0.349 notebook keeps the 1000 most-documented connectivities in the 5 ppm window; on held-out train
molecules that holds the answer for 43%. The 0.421 notebook adds a second pass chosen by the spectrum:
FPNet's f.z on the ECFP4 block alone (cheap: one Morgan fingerprint per structure) over the whole window,
the top N1 rows of that re-scored with the full fingerprint. This probe measures what that pass reaches,
before any notebook change or submission. v2 kept N1 = 5000; v3 keeps the whole window, so any N1 can be
read off locally.

Queries are a seeded sample of one fold-safe fold (`hold`, or the `np` natural products), scored by `fpnet_fold0.pt`, which never saw
them. Per query it writes every popularity-shortlist row and every pass-1 top-N1 row with its full f.z,
popularity, pass-1 rank and label, so recall and MRR@25 at any N1 / keep size are computed locally
afterwards. Writes `probe/part-NNNN.parquet` and `probe_stats.json` to /kaggle/working.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import random
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

N_QUERIES = 300  # capped at the fold's size
FOLD = "np"  # fold-safe split fold to sample: "hold" (v2, v3) or "np" (natural products, v4)
SEED = 1
N1: int | None = None  # pass-1 rows kept with full f.z; None = the whole window (v2 kept 5000, as the 0.421 notebook)
ECFP_SIZE = 4096  # raw fingerprint bits below this are the ECFP4 block
WORKERS = 4
CHECKPOINT = 50
OUT = Path("/kaggle/working")

W: dict[str, Any] = {}  # set before the pool forks; workers read it, never write it


def load_helpers() -> ModuleType:
    code = sorted(Path("/kaggle/input").glob("**/ranker_trainset.py"))
    if not code:
        raise SystemExit("ranker_trainset.py not found; is casmi-notebook-code attached?")
    sys.path.insert(0, str(code[0].parent))
    import ranker_trainset as rt

    return rt


def sample_fold(split: Path) -> list[str]:
    rows = pq.read_table(split, columns=["inchikey14", "fold"]).to_pylist()
    keys = sorted(r["inchikey14"] for r in rows if r["fold"] == FOLD)
    return random.Random(SEED).sample(keys, min(N_QUERIES, len(keys)))


def partial_scores(task: tuple[int, int, np.ndarray]) -> np.ndarray:
    """f.z on the ECFP4 block for tier rows [a, b); -inf where RDKit cannot parse the SMILES."""
    a, b, z_ecfp = task
    out = np.full(b - a, -np.inf, np.float32)
    for n, i in enumerate(range(a, b)):
        mol = W["Chem"].MolFromSmiles(W["tier"].smiles_at(i))
        if mol is not None:
            out[n] = W["gen"].GetFingerprintAsNumPy(mol)[W["ecfp_raw"]].astype(np.float32) @ z_ecfp
    return out


def packed_fingerprint(i: int) -> np.ndarray:
    fp = W["chem"].raw_fingerprint(W["tier"].smiles_at(i))
    bits = W["bits"]
    return np.packbits(fp[bits] if fp is not None else np.zeros(len(bits), np.uint8))


def spans(a: int, b: int, n: int) -> list[tuple[int, int]]:
    edges = np.linspace(a, b, n + 1).astype(int)
    return [(int(lo), int(hi)) for lo, hi in zip(edges[:-1], edges[1:]) if hi > lo]


def check_ecfp_block(fpnet: Any, tier: Any, a: int) -> None:
    """The ECFP4 block of the raw fingerprint must be the 4096-bit radius-2 Morgan fingerprint."""
    checked = 0
    for i in range(a, a + 20):
        smi = tier.smiles_at(i)
        raw, mol = fpnet.chem.raw_fingerprint(smi), W["Chem"].MolFromSmiles(smi)
        if raw is not None and mol is not None:
            assert np.array_equal(raw[:ECFP_SIZE], W["gen"].GetFingerprintAsNumPy(mol)), \
                "raw_fingerprint's first 4096 bits are not ECFP4: pass 1 would score the wrong bits"
            checked += 1
    assert checked, "no sampled tier row parsed; the ECFP4 block is unchecked"


def probe_rows(nb: ModuleType, pool: Any, key: str, mass: float, z: np.ndarray,
               stats: dict[str, Any]) -> list[dict[str, Any]]:
    """Shortlist and pass-1 rows of one query; [] when the window has no keyed row (counted as a miss)."""
    tier = W["tier"]
    tol = mass * nb.PPM * 1e-6
    a = int(np.searchsorted(tier.mass, mass - tol, "left"))
    b = int(np.searchsorted(tier.mass, mass + tol, "right"))
    keys = np.asarray(tier.ik14[a:b])
    stats["window_rows"].append(b - a)
    stats["truth_in_window"] += int(key.encode() in set(keys.tolist()))
    if not (keys != b"").any():
        return []
    pop = np.add(tier.lsid[a:b], tier.lpmid[a:b], dtype=np.float32)
    in_pop = set(nb.top_distinct(pop, keys, nb.SHORTLIST))
    z_ecfp = z[W["ecfp_cols"]].astype(np.float32)
    partial = np.concatenate(pool.map(partial_scores, [(lo, hi, z_ecfp) for lo, hi in spans(a, b, WORKERS * 4)]))
    order = np.argsort(-partial, kind="stable")[:N1]
    pass1_rank = {int(j): r for r, j in enumerate(order) if np.isfinite(partial[j])}
    picked = sorted((in_pop | set(pass1_rank)) & set(np.flatnonzero(keys != b"").tolist()))
    packed = pool.map(packed_fingerprint, [a + j for j in picked], chunksize=64)
    fz = np.unpackbits(np.stack(packed), axis=1)[:, :len(W["bits"])].astype(np.float32) @ z
    target = key.encode()
    return [dict(qid=key, key=keys[j].decode(), pop=float(pop[j]), fz=float(fz[n]),
                 pass1_rank=int(pass1_rank.get(j, -1)), in_pop=j in in_pop, label=bool(keys[j] == target))
            for n, j in enumerate(picked)]


def write_part(rows: list[dict[str, Any]], part: int, stats: dict[str, Any]) -> None:
    (OUT / "probe").mkdir(exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), OUT / "probe" / f"part-{part:04d}.parquet")
    (OUT / "probe_stats.json").write_text(json.dumps(stats))


def main() -> None:
    t0 = time.time()
    rt = load_helpers()
    nb = rt.load_notebook()
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    tier = nb.find_pubchem()
    fpnet = nb.FPNetScorer(rt.find("code/casmi/fpnet.py").parent.parent,
                           [rt.find("casmi26-fold-safe-fpnet/**/fpnet_fold0.pt")],
                           rt.find("casmi26-v2-pool/**/fp_bits.npy"))
    cols = np.flatnonzero(fpnet.bits < ECFP_SIZE)
    W.update(tier=tier, chem=fpnet.chem, bits=fpnet.bits, Chem=Chem, ecfp_cols=cols, ecfp_raw=fpnet.bits[cols],
             gen=rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=ECFP_SIZE))
    check_ecfp_block(fpnet, tier, len(tier.mass) // 2)
    keys = sample_fold(rt.find("casmi26-fold-safe-fpnet/**/split.parquet"))
    mols = rt.query_spectra(nb, rt.find("train.parquet"), set(keys))
    print(f"queries: {len(keys)} chosen, {len(mols)} with spectra; ECFP4 block checked  ({time.time()-t0:.0f}s)",
          flush=True)
    stats: dict[str, Any] = {"queries": 0, "no_mass": 0, "truth_in_window": 0, "window_rows": []}
    rows: list[dict[str, Any]] = []
    part = 0
    with mp.get_context("fork").Pool(WORKERS) as pool:
        for n, (key, mol) in enumerate(mols.items(), 1):
            mass = nb.neutral_mass(mol["mz"], mol["adduct"]) if mol["mz"] else None
            if mass is None:
                stats["no_mass"] += 1
                continue
            rows += probe_rows(nb, pool, key, mass, fpnet.logits(mol["records"]), stats)
            stats["queries"] += 1
            if n % 10 == 0:
                print(f"  {n}/{len(mols)}  in window {stats['truth_in_window']}/{stats['queries']}  "
                      f"rows {len(rows):,}  ({time.time()-t0:.0f}s)", flush=True)
            if n % CHECKPOINT == 0:
                write_part(rows, part, stats)
                rows, part = [], part + 1
    stats["seconds"] = round(time.time() - t0)
    if rows:
        write_part(rows, part, stats)
    (OUT / "probe_stats.json").write_text(json.dumps(stats))
    print(f"done: queries {stats['queries']}, truth in window {stats['truth_in_window']}, "
          f"median window {int(np.median(stats['window_rows']))} rows, {stats['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
