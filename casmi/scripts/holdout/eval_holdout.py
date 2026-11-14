"""Score the eval split with pool-style candidate lists, compare scorers and blends by MRR@25, and judge the gate.

Usage:
  eval_holdout.py lists                       build runs/lists.parquet (same-formula pool candidates)
  eval_holdout.py score <name> ...            score every list with runs/<name>/model via the kit's inference.Scorer
  eval_holdout.py compare <base> <cand> ...   paired MRR@25; a candidate "a+b" is a per-molecule z-blend
  eval_holdout.py gate                        the calibration gate: kit beats long, and kit beats kit+s1+s2+s3

Each list is the true structure plus up to 255 same-formula structures from pool_popularity.csv, most popular first;
the pipeline draws its candidates from the same pool. Ranks use pessimistic ties, as the kit's train.validate does.
`stats.py` is casmi/src/casmi/stats.py, shipped alongside.
"""

import json
import sys
from functools import cache
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from stats import paired

from common import ROOT, RUNS, TRAIN, WORK, load_kit, log

POOL = ROOT / "pool_popularity.csv"
LISTS = RUNS / "lists.parquet"
MAX_DECOYS, WORKERS = 255, 12
SPECTRUM_COLUMNS = {
    "ms2_mzs": "mzs",
    "ms2_normalized_intensities": "intensities",
    "precursor_mz": "precursor_mz",
    "adduct": "adduct",
    "ionization_mode": "ionization_mode",
    "collision_energy_ev": "collision_energy",
}


def featurise(smiles: list[str]) -> list[tuple[str, str, bytes]]:
    load_kit()
    from features import structure

    keyed: dict[str, tuple[str, bytes]] = {}
    for s in smiles:
        try:
            k, _, _, _, fp = structure(s)
        except (ValueError, RuntimeError):
            continue
        keyed.setdefault(k, (s, np.packbits(fp).tobytes()))
        if len(keyed) > MAX_DECOYS:
            break
    return [(k, s, fp) for k, (s, fp) in keyed.items()]


def build_lists() -> None:
    mol = pd.read_parquet(WORK / "molecules.parquet")
    fps = np.load(WORK / "fingerprints.npy", mmap_mode="r")
    ev = mol[mol.split.eq("eval")]
    pool = pd.read_csv(POOL, usecols=["smiles", "formula", "pop"])
    pool = pool[pool.formula.isin(set(ev.formula))].sort_values("pop", ascending=False)
    groups = {f: g.smiles.tolist() for f, g in pool.groupby("formula", sort=False)}
    with Pool(WORKERS) as p:
        decoys = dict(zip(groups, p.map(featurise, groups.values())))
    rows = []
    for m in ev.itertuples():
        found = [d for d in decoys.get(m.formula, []) if d[0] != m.key][:MAX_DECOYS]
        if found:
            rows += [(m.key, k, s, False, fp) for k, s, fp in found]
            rows.append((m.key, m.key, m.smiles, True, np.packbits(fps[m.Index]).tobytes()))
    out = pd.DataFrame(rows, columns=["mol_key", "cand_key", "cand_smiles", "is_true", "fp"])
    LISTS.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(LISTS, index=False)
    log(
        "LISTS",
        out.mol_key.nunique(),
        "molecules,",
        len(out),
        "rows; no same-formula decoys:",
        len(ev) - out.mol_key.nunique(),
    )


def eval_spectra(keys: list[str]) -> dict[str, list[dict[str, object]]]:
    """Raw train.parquet spectra per eval molecule, matched on inchikey14 (train.parquet has no spectrum ids)."""
    table = ds.dataset(TRAIN).to_table(
        columns=["inchikey14", *SPECTRUM_COLUMNS], filter=ds.field("inchikey14").isin(keys)
    )
    out: dict[str, list[dict[str, object]]] = {}
    for r in table.to_pandas().rename(columns=SPECTRUM_COLUMNS).to_dict("records"):
        out.setdefault(r.pop("inchikey14"), []).append(r)
    return out


def score(names: list[str]) -> None:
    load_kit()
    from inference import Scorer

    lists = pd.read_parquet(LISTS)
    spectra = eval_spectra(lists.mol_key.unique().tolist())
    lists = lists[lists.mol_key.isin(spectra)].reset_index(drop=True)
    fps = np.unpackbits(np.frombuffer(b"".join(lists.fp), np.uint8)).reshape(len(lists), -1)
    shared = dict(zip(lists.cand_smiles, fps))
    log("SPECTRA", len(spectra), "molecules matched on inchikey14,", sum(map(len, spectra.values())), "spectra")
    for name in names:
        scorer = Scorer(RUNS / name / "model", checkpoint="best_validation_model.pt")
        scorer.cache = shared
        scores = np.empty(len(lists), np.float32)
        for key, g in lists.groupby("mol_key", sort=False):
            scores[g.index.to_numpy()] = scorer.score_candidates(spectra[key], g.cand_smiles.tolist())
        lists[["mol_key", "cand_key", "is_true"]].assign(score=scores).to_parquet(RUNS / name / "eval_scores.parquet")
        log("SCORED", name, lists.mol_key.nunique(), "molecules")


@cache
def scores_of(name: str) -> pd.DataFrame:
    return pd.read_parquet(RUNS / name / "eval_scores.parquet")


@cache
def zscores_of(name: str) -> np.ndarray:
    f = scores_of(name)
    g = f.groupby("mol_key")["score"]
    return ((f.score - g.transform("mean")) / g.transform("std").replace(0, 1)).to_numpy()


def reciprocal_ranks(spec: str) -> pd.Series:
    names = spec.split("+")
    f = scores_of(names[0])[["mol_key", "is_true"]].assign(
        score=scores_of(names[0]).score if len(names) == 1 else np.mean([zscores_of(n) for n in names], axis=0)
    )
    true = f[f.is_true].set_index("mol_key").score
    better = f.assign(t=f.mol_key.map(true)).query("not is_true and score >= t").groupby("mol_key").size()
    rank = 1 + better.reindex(true.index, fill_value=0)
    return (1 / rank).where(rank <= 25, 0.0)


def compare(base: str, cands: list[str]) -> dict[str, dict[str, object]]:
    lib = pd.read_parquet(WORK / "molecules.parquet", columns=["key", "library"]).set_index("key").library
    b = reciprocal_ranks(base)
    out = {}
    for c in cands:
        r = reciprocal_ranks(c).reindex(b.index)
        row = paired(r.to_numpy(), b.to_numpy())
        row["by_library"] = {k: paired(r[g.index].to_numpy(), g.to_numpy())["delta"] for k, g in b.groupby(lib)}
        out[c] = row
        log(f"{c} vs {base}", json.dumps(row))
    return out


def gate() -> None:
    rows = compare("kit", ["long", "kit+s1+s2+s3"])
    checks = {"kit beats long": rows["long"]["delta"] < 0, "kit beats kit+s1+s2+s3": rows["kit+s1+s2+s3"]["delta"] < 0}
    verdict = dict(checks=checks, passed=all(checks.values()), rows=rows)
    (ROOT / "calibration_gate.json").write_text(json.dumps(verdict, indent=2))
    log("GATE", "PASS" if verdict["passed"] else "FAIL", json.dumps(checks))
    sys.exit(0 if verdict["passed"] else 1)


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "lists":
        build_lists()
    elif cmd == "score":
        score(args)
    elif cmd == "compare":
        compare(args[0], args[1:])
    else:
        gate()
