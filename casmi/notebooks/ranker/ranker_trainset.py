"""Kaggle kernel: the ranker's training set — features and labels for train molecules no model saw.

Queries are the fold-safe split's held-out molecules (casmi26-fold-safe-fpnet `split.parquet`): all
265 of the natural-product panel and a seeded sample of the `hold` fold. Their FPNet scores come from
`fpnet_fold0.pt`, trained without them, so a label is never something the scorer memorised.

Candidates are the production PubChem shortlist only. The train window is left out on purpose: a train
query's own structure always sits in it, so including it would teach the ranker to prefer train
structures, which at test time are almost never the answer. A query whose structure is not in the
shortlist has no positive; it is kept and counted, because that rate is the shortlist's recall.

Features are the ones the notebook can compute at test time for a shortlist candidate, all relative
within the molecule where the scale could shift: FPNet is the fold-0 checkpoint here and the two-checkpoint
ensemble at test time, so only its z-score and rank are features. `handset` is today's score, kept to
compare against, not to train on.

The notebook's code is imported from the `casmi-notebook-code` dataset rather than copied. Writes
`trainset/part-NNNN.parquet` (one per CHECKPOINT queries) and `trainset_stats.json` to /kaggle/working.
"""

from __future__ import annotations

import importlib.util
import json
import random
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

N_HOLD = 1735
SEED = 0
CHECKPOINT = 250
FAIL_FAST_AFTER = 50  # queries; all-zero labels by then means the keys do not match
INPUT = Path("/kaggle/input")
OUT = Path("/kaggle/working")


def find(pattern: str) -> Path:
    hits = sorted(INPUT.glob(f"**/{pattern}"))
    if not hits:
        raise SystemExit(f"{pattern} not found under {INPUT}; is its dataset attached?")
    return hits[0]


def load_notebook():
    spec = importlib.util.spec_from_file_location("nb", find("kaggle_submission.py"))
    nb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(nb)  # installs rdkit from the attached wheel, as the notebook does
    return nb


def queries(split: Path) -> dict[str, str]:
    """inchikey14 -> fold for every np-panel molecule and a seeded sample of the hold fold."""
    t = pq.read_table(split, columns=["inchikey14", "fold"]).to_pylist()
    folds = {r["inchikey14"]: r["fold"] for r in t}
    hold = sorted(k for k, f in folds.items() if f == "hold")
    picked = random.Random(SEED).sample(hold, N_HOLD)
    out = {k: folds[k] for k in picked} | {k: f for k, f in folds.items() if f == "np"}
    assert all(len(k) == 14 and k.isupper() for k in out), "split inchikey14 is not an InChIKey block"
    return out


def query_spectra(nb, train: Path, keys: set[str]) -> dict[str, dict]:
    wanted = pa.array(sorted(keys))
    rows = []
    for batch in pq.ParquetFile(train).iter_batches(batch_size=200_000, columns=["inchikey14", *nb.SPECTRUM_COLUMNS]):
        rows += batch.filter(pc.is_in(batch.column("inchikey14"), value_set=wanted)).to_pylist()
    return nb.group_spectra(rows, "inchikey14")


class Prior:
    """The popularity prior's extra arrays, row-aligned with the PubChem tier."""

    def __init__(self, root: Path) -> None:
        self.lpatent = np.load(root / "pc_lpatent.npy", mmap_mode="r")
        self.np_bits = np.load(root / "pc_np.npy", mmap_mode="r")


def candidate_rows(nb, tier, prior: Prior, key: str, mol: dict, mass: float, fpnet) -> list[dict]:
    tol = mass * nb.PPM * 1e-6
    a = int(np.searchsorted(tier.mass, mass - tol, "left"))
    b = int(np.searchsorted(tier.mass, mass + tol, "right"))
    keys = np.asarray(tier.ik14[a:b])
    pop = np.add(tier.lsid[a:b], tier.lpmid[a:b], dtype=np.float32)
    picks = nb.top_distinct(pop, keys, nb.SHORTLIST)
    if not picks:
        return []
    rows = [a + j for j in picks]
    smiles = [tier.smiles_at(i) for i in rows]
    fz = fpnet.scores(mol["records"], smiles)
    fpnet.forget()  # candidates rarely repeat across queries; the cache would only grow
    pops = pop[picks]
    formulas = [nb.ForwardRescorer.formula(s) for s in smiles]
    counts = {f: formulas.count(f) for f in set(formulas)}
    handset = nb.blend(fz, pops)
    fz_rank = np.argsort(np.argsort(-fz))
    pop_rank = np.argsort(np.argsort(-pops))
    sd = fz.std() or 1.0
    return [dict(
        qid=key, label=int(keys[j].decode() == key),
        fz_z=float((fz[n] - fz.mean()) / sd), fz_rank=int(fz_rank[n]),
        lsid=float(tier.lsid[i]), lpmid=float(tier.lpmid[i]), lpatent=float(prior.lpatent[i]),
        pop=float(pops[n]), pop_rank=int(pop_rank[n]),
        np_coconut=int(prior.np_bits[i] & 1 > 0), np_lotus=int(prior.np_bits[i] & 2 > 0),
        np_npatlas=int(prior.np_bits[i] & 4 > 0),
        ppm_error=float((tier.mass[i] - mass) / mass * 1e6),
        n_isomers=counts[formulas[n]], n_candidates=len(rows), handset=float(handset[n]))
        for n, (i, j) in enumerate(zip(rows, picks))]


def record(stats: dict, fold: str, hit: bool) -> None:
    fs = stats["by_fold"].setdefault(fold, {"queries": 0, "with_positive": 0})
    for d in (stats, fs):
        d["queries"] += 1
        d["with_positive"] += int(hit)


def write_part(rows: list[dict], part: int, stats: dict) -> None:
    (OUT / "trainset").mkdir(exist_ok=True)
    table = pa.Table.from_pylist(rows)
    table = table.cast(pa.schema([pa.field(f.name, pa.float32() if pa.types.is_floating(f.type) else
                                           pa.int16() if pa.types.is_integer(f.type) else f.type)
                                  for f in table.schema]))
    pq.write_table(table, OUT / "trainset" / f"part-{part:04d}.parquet")
    (OUT / "trainset_stats.json").write_text(json.dumps(stats, indent=1))


def main() -> None:
    t0 = time.time()
    nb = load_notebook()
    tier = nb.find_pubchem()
    prior = Prior(find("pc_lpatent.npy").parent)
    fpnet = nb.FPNetScorer(find("code/casmi/fpnet.py").parent.parent,
                           [find("casmi26-fold-safe-fpnet/**/fpnet_fold0.pt")],
                           find("casmi26-v2-pool/**/fp_bits.npy"))
    folds = queries(find("casmi26-fold-safe-fpnet/**/split.parquet"))
    mols = query_spectra(nb, find("train.parquet"), set(folds))
    print(f"queries: {len(folds)} chosen, {len(mols)} with spectra  ({time.time()-t0:.0f}s)", flush=True)
    rows: list[dict] = []
    part = 0
    stats = {"queries": 0, "with_positive": 0, "no_mass": 0, "by_fold": {}}
    for n, (key, mol) in enumerate(mols.items(), 1):
        mass = nb.neutral_mass(mol["mz"], mol["adduct"]) if mol["mz"] else None
        if mass is None:
            stats["no_mass"] += 1
        else:
            cand = candidate_rows(nb, tier, prior, key, mol, mass, fpnet)
            rows += [{**r, "fold": folds[key]} for r in cand]
            record(stats, folds[key], any(r["label"] for r in cand))
        if n % 25 == 0:
            print(f"  {n}/{len(mols)}  recall {stats['with_positive']}/{stats['queries']}  "
                  f"rows {len(rows):,}  ({time.time()-t0:.0f}s)", flush=True)
        if stats["queries"] >= FAIL_FAST_AFTER and not stats["with_positive"]:
            raise SystemExit("no query has its own structure among the candidates: the keys do not match")
        if n % CHECKPOINT == 0 and rows:
            write_part(rows, part, stats)
            rows, part = [], part + 1
    stats["seconds"] = round(time.time() - t0)
    if rows:
        write_part(rows, part, stats)
    (OUT / "trainset_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"done: {json.dumps(stats)}", flush=True)


if __name__ == "__main__":
    main()
