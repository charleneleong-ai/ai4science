"""Paired comparison of two rankers on the same molecules."""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

BOOTSTRAP = 2000


def paired(model_rr: npt.NDArray[np.float64], base_rr: npt.NDArray[np.float64], seed: int = 0) -> dict[str, Any]:
    """Per-molecule reciprocal ranks of a change and its baseline: MRR of both, the mean paired difference with a bootstrap 95% interval, and win / loss counts."""
    diff = model_rr - base_rr
    rng = np.random.default_rng(seed)
    boot = diff[rng.integers(0, len(diff), (BOOTSTRAP, len(diff)))].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"queries": len(diff), "base": float(base_rr.mean()), "model": float(model_rr.mean()),
            "delta": float(diff.mean()), "ci95": [float(lo), float(hi)],
            "better": int((diff > 0).sum()), "worse": int((diff < 0).sum())}
