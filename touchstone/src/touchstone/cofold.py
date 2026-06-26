"""Co-fold cross-check: corroborate a design against an INDEPENDENT predictor.

The geometry and physics tiers judge a design's own coordinates. This stage asks a
different question: if a second, unrelated structure predictor (Chai-1, AllMetal3D,
…) folds the same sequence with the same metal, does it place the metal in the
coordination environment the design claims? Two predictors agreeing is the co-fold
analogue of two MLIPs agreeing — independent corroboration, not self-consistency.

Predictor-agnostic: the verifier consumes a predicted CoordinationSite via a
`provider` callback, not any particular model. The heavy inference (running Chai-1
on a GPU box) lives in scripts/chai_crosscheck.py and writes a structure the
provider parses.
"""

from __future__ import annotations

from collections import Counter
from typing import Callable

import numpy as np

from .core import BinderDesign, CoordinationSite, Verdict
from .geometry.parse import coordination_site


def _donor_overlap(a: tuple[str, ...], b: tuple[str, ...]) -> float:
    """Fraction of the larger donor multiset matched element-for-element."""
    if not a or not b:
        return 0.0
    return sum((Counter(a) & Counter(b)).values()) / max(len(a), len(b))


def cofold_agreement(
    reference: CoordinationSite,
    predicted: CoordinationSite,
    cn_tol: int = 1,
    bond_tol: float = 0.4,
    donor_tol: float = 0.5,
) -> Verdict:
    """Score how well an independent prediction recovers the design's coordination.

    Trusts when the predicted site matches on coordination number, donor identity,
    and mean bond length. Strong disagreement is `weak` (not corroborated, but
    judgeable); a predictor that places no metal site at all is `defer`.
    """
    if predicted.is_empty:
        return Verdict.defer("co-fold placed no coordinating atoms")

    cn_delta = abs(reference.coordination_number - predicted.coordination_number)
    overlap = _donor_overlap(reference.ligand_elems, predicted.ligand_elems)
    bond_delta = abs(float(reference.bond_lengths().mean()) - float(predicted.bond_lengths().mean()))

    score = float(overlap * np.exp(-0.5 * (cn_delta / cn_tol) ** 2) * np.exp(-0.5 * (bond_delta / bond_tol) ** 2))
    agree = cn_delta <= cn_tol and overlap >= donor_tol and bond_delta <= bond_tol
    reason = (
        f"co-fold CN {predicted.coordination_number} vs {reference.coordination_number}, "
        f"donor overlap {overlap:.0%}, Δbond {bond_delta:.2f} Å"
    )
    return Verdict(score, trust=agree, ood=False, reason=reason)


class CofoldCrossCheck:
    """Verifier over an independent co-fold prediction. `provider` maps a design to
    its predicted CoordinationSite (or None when no prediction is available)."""

    def __init__(
        self,
        provider: Callable[[BinderDesign], CoordinationSite | None],
        *,
        cn_tol: int = 1,
        bond_tol: float = 0.4,
        donor_tol: float = 0.5,
    ):
        self.provider = provider
        self.cn_tol = cn_tol
        self.bond_tol = bond_tol
        self.donor_tol = donor_tol

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            predicted = self.provider(design)
        except Exception as e:  # inference / parse failure ⇒ can't cross-check
            return Verdict.defer(f"co-fold prediction failed: {type(e).__name__}")
        if predicted is None:
            return Verdict.defer("no co-fold prediction available")
        return cofold_agreement(design.site, predicted, self.cn_tol, self.bond_tol, self.donor_tol)


def cif_provider(predictions: dict[str, str], metal_atom: str = "NI", metal: str = "Ni2+"):
    """File-based provider: parse each design's predicted CIF/PDB from a
    {design.source: predicted_path} map (e.g. Chai-1 / AllMetal3D outputs)."""

    def provide(design: BinderDesign) -> CoordinationSite | None:
        path = predictions.get(design.source or "")
        return coordination_site(path, metal_atom, metal) if path else None

    return provide
