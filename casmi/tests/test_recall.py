from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from casmi.recall import Config, Query, load, reciprocal_ranks


@pytest.fixture
def query() -> Query:
    """Shortlist {A, B}; pass 1 also found C twice (stereo rows) and the answer D, which is unpopular."""
    return Query(key=np.array(["A", "B", "C", "C", "D"], dtype=object),
                 pop=np.array([9, 8, 0, 0, 0], np.float32), fz=np.array([0, 1, 5, 4, 3], np.float32),
                 pass1_rank=np.array([-1, 3, 0, 1, 2], np.int32),
                 in_pop=np.array([1, 1, 0, 0, 0], bool), label=np.array([0, 0, 0, 0, 1], bool))


class TestSecondPass:
    @pytest.mark.parametrize(("cfg", "keys"), [
        (Config(0, 0), ["A", "B"]),
        (Config(5000, 1), ["A", "B", "C"]),  # one row per connectivity, the best f.z
        (Config(5000, 2), ["A", "B", "C", "D"]),
        (Config(2, 5), ["A", "B", "C"]),  # D's pass-1 rank 2 is beyond N1 = 2
    ])
    def test_candidates(self, query: Query, cfg: Config, keys: list[str]) -> None:
        rows = query.candidates(cfg)
        assert query.key[rows].tolist() == keys and query.fz[rows[2:]].tolist() == sorted(query.fz[rows[2:]], reverse=True)

    def test_answer_reached_only_by_the_second_pass_is_ranked(self, query: Query) -> None:
        assert query.reciprocal_rank(Config(0, 0)) == (False, 0.0)
        hit, rr = query.reciprocal_rank(Config(5000, 2))
        assert hit and 0 < rr < 1  # popularity keeps A and B above it

    def test_queries_without_rows_count_as_misses(self, tmp_path: Path, query: Query) -> None:
        rows = [dict(qid="Q", key=str(k), pop=float(p), fz=float(f), pass1_rank=int(r), in_pop=bool(i), label=bool(lb))
                for k, p, f, r, i, lb in zip(query.key, query.pop, query.fz, query.pass1_rank, query.in_pop, query.label)]
        pq.write_table(pa.Table.from_pylist(rows), tmp_path / "part-0000.parquet")
        hits, rr = reciprocal_ranks(load(tmp_path), Config(5000, 2), total=3)
        assert hits == 1 and len(rr) == 3 and rr[1:].tolist() == [0.0, 0.0]
