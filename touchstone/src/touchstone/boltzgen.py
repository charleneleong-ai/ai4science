"""Read BoltzGen fold-output confidence (`.npz`) — the generator grading its own prediction.

Shared by the RLVR tooling (`scripts/rlvr_select.py`) and the confidence-vs-verdict view
(`scripts/boltzgen_scores.py`): both pull the same iPTM/pLDDT/pTM from a design's fold-output
`.npz`, taking the best refold sample.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def best(npz, key: str) -> float | None:
    """Max of a per-refold-sample confidence array, or None if the key is absent."""
    if key not in npz:
        return None
    return float(np.atleast_1d(npz[key]).astype(float).max())


def boltzgen_confidence(npz_path: str | Path) -> dict | None:
    """`{iptm, plddt, ptm}` for a design's BoltzGen fold-output `.npz` (best refold sample),
    or None if the file is absent. iPTM falls back across BoltzGen's key names."""
    npz_path = Path(npz_path)
    if not npz_path.exists():
        return None
    npz = np.load(npz_path, allow_pickle=True)
    return {
        "iptm": best(npz, "design_to_target_iptm") or best(npz, "ligand_iptm") or best(npz, "iptm"),
        "plddt": best(npz, "complex_plddt"),
        "ptm": best(npz, "design_ptm") or best(npz, "ptm"),
    }
