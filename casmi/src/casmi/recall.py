"""What a spectrum-driven second PubChem pass adds, read from the recall-probe kernel's rows.

Per query the probe kept every popularity-shortlist row and every pass-1 top-N1 row with its full
f.z. A configuration (N1, keep) adds the `keep` best distinct connectivities by full f.z among the
pass-1 top N1 to the shortlist; candidates are then ranked by the notebook's blend (z(f.z) + 0.25 ×
popularity, z over the molecule's candidates). Ties count against the answer, as in `ranker`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import pyarrow.parquet as pq

POP_WEIGHT = 0.25  # the notebook's blend
TOP_K = 25


@dataclass(frozen=True)
class Config:
    n1: int  # pass-1 rows re-scored with the full fingerprint; 0 = popularity shortlist only
    keep: int  # distinct connectivities taken from them


@dataclass
class Query:
    key: npt.NDArray[np.object_]
    pop: npt.NDArray[np.float32]
    fz: npt.NDArray[np.float32]
    pass1_rank: npt.NDArray[np.int32]
    in_pop: npt.NDArray[np.bool_]
    label: npt.NDArray[np.bool_]

    def candidates(self, cfg: Config) -> npt.NDArray[np.int64]:
        """Row indices: the shortlist plus the `keep` best new connectivities by f.z, one row each."""
        chosen = list(np.flatnonzero(self.in_pop))
        seen = set(self.key[chosen])
        if cfg.keep:
            pool = np.flatnonzero((self.pass1_rank >= 0) & (self.pass1_rank < cfg.n1))
            added = 0
            for i in pool[np.argsort(-self.fz[pool], kind="stable")]:
                if self.key[i] in seen:
                    continue
                seen.add(self.key[i])
                chosen.append(i)
                added += 1
                if added == cfg.keep:
                    break
        return np.array(chosen, dtype=np.int64)

    def reciprocal_rank(self, cfg: Config) -> tuple[bool, float]:
        """Whether the answer is a candidate, and 1 / its rank under the blend (0 beyond the top 25)."""
        rows = self.candidates(cfg)
        hit = np.flatnonzero(self.label[rows])
        if not hit.size:
            return False, 0.0
        fz = self.fz[rows]
        sd = fz.std() or 1.0
        score = (fz - fz.mean()) / sd + POP_WEIGHT * self.pop[rows]
        rank = int((score >= score[hit[0]]).sum())
        return True, 1.0 / rank if rank <= TOP_K else 0.0


def load(root: Path) -> list[Query]:
    parts = sorted(root.glob("part-*.parquet"))
    if not parts:
        raise FileNotFoundError(f"no part-*.parquet under {root}; has the probe kernel written any?")
    table = pa.concat_tables(pq.read_table(p) for p in parts).sort_by("qid")
    qid = table.column("qid").to_numpy(zero_copy_only=False)
    cols = {c: table.column(c).to_numpy(zero_copy_only=False) for c in
            ("key", "pop", "fz", "pass1_rank", "in_pop", "label")}
    edges = np.r_[np.flatnonzero(np.r_[True, qid[1:] != qid[:-1]]), len(qid)]
    return [Query(**{c: v[a:b] for c, v in cols.items()}) for a, b in zip(edges[:-1], edges[1:])]


def reciprocal_ranks(queries: list[Query], cfg: Config, total: int) -> tuple[int, npt.NDArray[np.float64]]:
    """Answers reached, and one reciprocal rank per query; queries with no rows at all count as 0."""
    results = [q.reciprocal_rank(cfg) for q in queries]
    rr = np.zeros(total)
    rr[:len(results)] = [r for _, r in results]
    return sum(h for h, _ in results), rr
