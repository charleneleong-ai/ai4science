"""LightGBM LambdaRank over the ranker training set, judged against the hand-set score.

Evaluation is grouped K-fold over query molecules: every candidate of a molecule lands in the same
fold, so no molecule is ranked by a model that saw it. Within each round, the next fold is the
early-stopping set and the rest train. The comparison is paired per molecule (same candidate lists,
same positives), with a bootstrap interval on the mean difference. It is a go / no-go for submitting,
not a leaderboard estimate: these are train molecules scored by the fold-0 FPNet.

With at most one positive per molecule, LightGBM's MAP@25 equals MRR@25 up to how it scores molecules
without a positive (a constant shift), so it serves as the early-stopping metric.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import lightgbm as lgb
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import pyarrow.parquet as pq

FEATURES = ["fz_z", "fz_rank", "lsid", "lpmid", "lpatent", "pop", "pop_rank", "np_coconut", "np_lotus",
            "np_npatlas", "ppm_error", "n_isomers", "n_candidates"]
TOP_K = 25
PARAMS: dict[str, Any] = {
    "objective": "lambdarank", "metric": "map", "eval_at": [TOP_K], "lambdarank_truncation_level": 40,
    "learning_rate": 0.05, "num_leaves": 31, "min_data_in_leaf": 200, "feature_fraction": 0.9,
    "bagging_fraction": 0.8, "bagging_freq": 1, "verbose": -1, "seed": 0,
}
BOOTSTRAP = 2000


@dataclass
class TrainSet:
    """Candidates sorted by molecule; `starts` are the offsets where each molecule begins."""

    X: npt.NDArray[np.float32]
    y: npt.NDArray[np.int32]
    handset: npt.NDArray[np.float32]
    qids: npt.NDArray[np.object_]  # one per molecule
    panel: npt.NDArray[np.object_]  # the split.parquet fold (hold / np), one per molecule
    starts: npt.NDArray[np.int64]

    @classmethod
    def load(cls, root: Path) -> TrainSet:
        table = pa.concat_tables(pq.read_table(p) for p in sorted(root.glob("part-*.parquet")))
        table = table.sort_by("qid")
        qid = table.column("qid").to_numpy(zero_copy_only=False)
        starts = np.flatnonzero(np.r_[True, qid[1:] != qid[:-1]])
        return cls(X=np.column_stack([table.column(f).to_numpy().astype(np.float32) for f in FEATURES]),
                   y=table.column("label").to_numpy().astype(np.int32),
                   handset=table.column("handset").to_numpy().astype(np.float32),
                   qids=qid[starts], panel=table.column("fold").to_numpy(zero_copy_only=False)[starts],
                   starts=starts)

    @property
    def sizes(self) -> npt.NDArray[np.int64]:
        return np.diff(np.r_[self.starts, len(self.y)])

    def rows(self, queries: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        """Candidate row indices of the given molecules, in molecule order."""
        return np.concatenate([np.arange(s, s + n) for s, n in zip(self.starts[queries], self.sizes[queries])])


def reciprocal_ranks(scores: npt.NDArray[np.floating], labels: npt.NDArray[np.integer],
                     starts: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    """Per molecule: 1 / rank of its positive if within the top 25, else 0. Ties count against it."""
    out = np.zeros(len(starts))
    for q, (a, b) in enumerate(zip(starts, np.r_[starts[1:], len(scores)])):
        pos = np.flatnonzero(labels[a:b])
        if pos.size:
            rank = int((scores[a:b] >= scores[a + pos[0]]).sum())
            out[q] = 1.0 / rank if rank <= TOP_K else 0.0
    return out


def query_folds(qids: npt.NDArray[np.object_], k: int) -> npt.NDArray[np.int64]:
    """A stable fold per molecule, from a hash of its key (independent of the sample drawn)."""
    return np.array([int(hashlib.md5(q.encode()).hexdigest(), 16) % k for q in qids])


def dataset(ts: TrainSet, queries: npt.NDArray[np.int64]) -> lgb.Dataset:
    rows = ts.rows(queries)
    return lgb.Dataset(ts.X[rows], ts.y[rows], group=ts.sizes[queries], feature_name=FEATURES,
                       free_raw_data=False)


def split(fold: npt.NDArray[np.int64], i: int, k: int) -> tuple[npt.NDArray[np.int64], ...]:
    """Round i's molecules: train (every other fold), early-stopping (fold i+1) and test (fold i)."""
    if k < 3:
        raise ValueError(f"need at least 3 folds for train / early-stopping / test, got {k}")
    stop = (i + 1) % k
    return np.flatnonzero((fold != i) & (fold != stop)), np.flatnonzero(fold == stop), np.flatnonzero(fold == i)


def cross_validate(ts: TrainSet, k: int, rounds: int) -> tuple[npt.NDArray[np.float32], list[int]]:
    """Out-of-fold ranker scores for every candidate, and each round's best iteration."""
    fold = query_folds(ts.qids, k)
    scores = np.zeros(len(ts.y), dtype=np.float32)
    best = []
    for i in range(k):
        train, stop, test = split(fold, i, k)
        model = lgb.train(PARAMS, dataset(ts, train), rounds, valid_sets=[dataset(ts, stop)],
                          callbacks=[lgb.early_stopping(100, verbose=False)])
        rows = ts.rows(test)
        scores[rows] = model.predict(ts.X[rows], num_iteration=model.best_iteration)
        best.append(model.best_iteration)
    return scores, best


def paired(model_rr: npt.NDArray[np.float64], base_rr: npt.NDArray[np.float64], seed: int = 0) -> dict[str, Any]:
    """MRR of both, the mean paired difference with a bootstrap 95% interval, and win / loss counts."""
    diff = model_rr - base_rr
    rng = np.random.default_rng(seed)
    boot = diff[rng.integers(0, len(diff), (BOOTSTRAP, len(diff)))].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"queries": len(diff), "handset": float(base_rr.mean()), "ranker": float(model_rr.mean()),
            "delta": float(diff.mean()), "ci95": [float(lo), float(hi)],
            "better": int((diff > 0).sum()), "worse": int((diff < 0).sum())}


def evaluate(ts: TrainSet, k: int, rounds: int) -> dict[str, Any]:
    scores, best = cross_validate(ts, k, rounds)
    model_rr = reciprocal_ranks(scores, ts.y, ts.starts)
    base_rr = reciprocal_ranks(ts.handset, ts.y, ts.starts)
    groups = {"all": np.ones(len(ts.qids), bool)} | {p: ts.panel == p for p in sorted(set(ts.panel))}
    return {"best_iterations": best,
            "groups": {g: paired(model_rr[m], base_rr[m]) for g, m in groups.items()}}


def fit_final(ts: TrainSet, rounds: int) -> lgb.Booster:
    return lgb.train(PARAMS, dataset(ts, np.arange(len(ts.qids))), rounds)


def save(model: lgb.Booster, report: dict[str, Any], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    model.save_model(str(out / "ranker.txt"))
    gain = model.feature_importance("gain")
    report["importance_gain"] = dict(sorted(zip(FEATURES, (gain / (gain.sum() or 1.0)).round(4).tolist()),
                                            key=lambda kv: -kv[1]))
    (out / "metrics.json").write_text(json.dumps(report, indent=1))
