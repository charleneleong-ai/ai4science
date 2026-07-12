"""Topology-reorganization on binding — is the metal site preorganized?

The MLIP dynamics verifier asks whether the *bound* site survives; this asks the complementary
design-quality question: does the fold already hold the coordination shell in place, or does the
shell only form around the metal? A preorganized site pays little entropic cost on binding
(stronger, more specific); a site that collapses without the metal is a weaker design.

A cheap MLIP proxy for the full apo/holo fold diff: remove the metal, relax the freed first-shell
donors against the frozen backbone (`relax_apo`), and measure how far they drift from the holo
positions. Reuses the same backbone + cluster plumbing as the static / dynamics / selectivity tiers.
"""

from __future__ import annotations

import numpy as np

from ..core import BinderDesign, Verdict
from .mlip import MLIPBase, relax_apo


class TrsVerifier(MLIPBase):
    """Trusts a preorganized site — first-shell donors barely move (≤ `trust_drift` Å) when the
    metal is removed and relaxed; `weak` when they reorganize moderately; `defer` when the shell
    fully collapses (> `ood_drift`) or the relaxation fails."""

    def __init__(self, backbone: str = "mace_mp", *, trust_drift: float = 1.0, ood_drift: float = 3.0,
                 steps: int = 200, **kw):
        super().__init__(backbone, **kw)
        self.trust_drift = trust_drift
        self.ood_drift = ood_drift
        self.steps = steps

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            r = relax_apo(self._cluster(design), self.calc, self._metal(design), self.cutoff, steps=self.steps)
        except Exception as e:  # cluster / relaxation failure ⇒ can't judge preorganization
            return Verdict.defer(f"apo relaxation failed: {type(e).__name__}")

        score = float(np.exp(-r.donor_drift))  # 1 at zero drift, decaying with reorganization
        metrics = {"donor_drift_angstrom": round(r.donor_drift, 3), "donors": r.donors}
        reason = f"donors drift {r.donor_drift:.2f} Å with the metal removed ({r.donors} donors)"
        if r.donor_drift > self.ood_drift:  # shell fully reorganizes ⇒ off the preorganized manifold
            return Verdict.defer(reason, score=score, metrics=metrics)
        return Verdict(score, trust=r.donor_drift <= self.trust_drift, ood=False, reason=reason, metrics=metrics)
