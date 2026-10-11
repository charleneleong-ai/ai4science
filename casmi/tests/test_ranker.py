from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ranker = pytest.importorskip("casmi.ranker")


def write_trainset(root: Path, n_queries: int = 60, n_cand: int = 40, seed: int = 0) -> None:
    """Molecules whose positive is the candidate with the highest `pop`; `handset` is noise."""
    rng = np.random.default_rng(seed)
    rows = []
    for q in range(n_queries):
        pos = rng.integers(n_cand)
        pop = rng.normal(size=n_cand)
        pop[pos] = pop.max() + 1
        for c in range(n_cand):
            rows.append({f: float(rng.normal()) for f in ranker.FEATURES} | {
                "qid": f"Q{q:013d}", "fold": "np" if q % 4 == 0 else "hold", "pop": float(pop[c]),
                "label": int(c == pos), "handset": float(rng.normal())})
    rng.shuffle(rows)  # load must regroup candidates by molecule
    pq.write_table(pa.Table.from_pylist(rows), root / "part-0000.parquet")


class TestReciprocalRanks:
    @pytest.mark.parametrize(("scores", "labels", "expected"), [
        ([0.9, 0.5, 0.1], [0, 1, 0], 0.5),
        ([0.5, 0.5, 0.1], [0, 1, 0], 0.5),  # a tie counts against the positive
        ([0.1, 0.2, 0.3], [0, 0, 0], 0.0),  # no positive in the shortlist
    ])
    def test_single_molecule(self, scores: list[float], labels: list[int], expected: float) -> None:
        rr = ranker.reciprocal_ranks(np.array(scores), np.array(labels), np.array([0]))
        assert rr.tolist() == [expected]

    def test_beyond_top_k_scores_zero(self) -> None:
        scores = np.arange(30, 0, -1, dtype=float)
        labels = np.zeros(30, int)
        labels[25] = 1
        assert ranker.reciprocal_ranks(scores, labels, np.array([0])).tolist() == [0.0]


class TestTrainSet:
    def test_load_groups_candidates_by_molecule(self, tmp_path: Path) -> None:
        write_trainset(tmp_path, n_queries=5, n_cand=7)
        ts = ranker.TrainSet.load(tmp_path)
        assert len(ts.qids) == 5 and ts.sizes.tolist() == [7] * 5
        assert all(ts.y[ts.rows(np.array([q]))].sum() == 1 for q in range(5))

    def test_folds_ignore_input_order(self) -> None:
        qids = np.array([f"Q{i:013d}" for i in range(200)])
        folds = ranker.query_folds(qids, 5)
        assert np.array_equal(folds, ranker.query_folds(qids[::-1], 5)[::-1])
        assert set(folds) == set(range(5))

    @pytest.mark.parametrize("i", range(4))
    def test_split_is_disjoint_and_complete(self, i: int) -> None:
        fold = np.arange(40) % 4
        train, stop, test = ranker.split(fold, i, 4)
        assert set(fold[test]) == {i} and set(fold[stop]) == {(i + 1) % 4}
        assert sorted(np.r_[train, stop, test]) == list(range(40))

    def test_split_needs_three_folds(self) -> None:
        with pytest.raises(ValueError, match="at least 3"):
            ranker.split(np.arange(10) % 2, 0, 2)


def test_ranker_learns_the_signal_the_baseline_misses(tmp_path: Path) -> None:
    write_trainset(tmp_path)
    report = ranker.evaluate(ranker.TrainSet.load(tmp_path), k=3, rounds=50)
    overall = report["groups"]["all"]
    assert overall["ranker"] > 0.9 > 0.3 > overall["handset"]
    assert overall["ci95"][0] > 0 and set(report["groups"]) == {"all", "hold", "np"}
