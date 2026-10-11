"""Read BoltzGen fold-output confidence (`.npz`) — the generator grading its own prediction.

Shared by the RLVR tooling (`scripts/rlvr_select.py`) and the confidence-vs-verdict view
(`scripts/boltzgen_scores.py`): both pull the same iPTM/pLDDT/pTM from a design's fold-output
`.npz`, which stores one value per refold (diffusion) sample.

The reduction over those samples defaults to the **mean**, not the best. Measured over 57
designs folded with 3 samples each (tnf-binder, PR #81), best-of-3 and the mean differ by up
to 0.35 *within a single design*, and the induced ranking reorders substantially — a design
ranked 1st on best-of-3 fell mid-pack on the mean. Ranking on the best sample selects for
sampling variance, not design quality: n samples give the generator n chances to flatter one
design. `iptm_best` rides along next to the chosen aggregate so the gap stays visible.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Literal

import numpy as np

Aggregate = Literal["mean", "best"]

REDUCERS: dict[Aggregate, Callable[[np.ndarray], np.floating]] = {"mean": np.mean, "best": np.max}

# BoltzGen renamed its confidence outputs across versions, so iPTM and pTM are fallback chains
# over *key names* — the first key the file actually carries wins.
IPTM_KEYS = ("design_to_target_iptm", "ligand_iptm", "iptm")
PTM_KEYS = ("design_ptm", "ptm")


def samples(npz, keys: tuple[str, ...]) -> np.ndarray | None:
    """The per-refold-sample array under the first of `keys` the file carries, or None if it
    carries none of them. Keyed on presence, not truthiness — a real 0.0 iPTM is a verdict,
    not a missing key."""
    for key in keys:
        if key in npz:
            return np.atleast_1d(np.asarray(npz[key], dtype=float))
    return None


def reduce_samples(values: np.ndarray | None, aggregate: Aggregate = "mean") -> float | None:
    """One number from a per-sample array, passing None through for an absent key."""
    return None if values is None else float(REDUCERS[aggregate](values))


def boltzgen_confidence(
    npz_path: str | Path, aggregate: Aggregate = "mean"
) -> dict[str, float | None] | None:
    """`{iptm, plddt, ptm}` for a design's BoltzGen fold-output `.npz`, each reduced over the
    refold samples by `aggregate` ("mean" by default, or "best"), plus `iptm_best` — so a caller
    ranking on the mean can still see how much of the headline number was sampling luck. None if
    the file is absent."""
    npz_path = Path(npz_path)
    if not npz_path.exists():
        return None
    npz = np.load(npz_path, allow_pickle=True)
    iptm = samples(npz, IPTM_KEYS)  # read once, reduced twice
    return {
        "iptm": reduce_samples(iptm, aggregate),
        "plddt": reduce_samples(samples(npz, ("complex_plddt",)), aggregate),
        "ptm": reduce_samples(samples(npz, PTM_KEYS), aggregate),
        "iptm_best": reduce_samples(iptm, "best"),
    }
