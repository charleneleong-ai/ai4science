"""The asset: a generator-blind verifier of metal-coordination geometry.

Scores a design's coordination site against the reference distribution for its metal,
and flags sites that sit off the manifold (implausible, or pushed there by extreme
conditions). Independent of whichever generator produced the design.
"""

from __future__ import annotations

import numpy as np

from ..core import BinderDesign, Verdict
from .reference import MockReference, ReferenceDistribution


class GeometryVerifier:
    def __init__(
        self,
        reference: ReferenceDistribution | None = None,
        trust_z: float = 2.0,
        ood_z: float = 3.0,
    ):
        self.reference = reference or MockReference()
        self.trust_z = trust_z  # bond-length z within this ⇒ geometry trusted
        self.ood_z = ood_z  # bond-length z beyond this ⇒ off-manifold ⇒ defer

    def verify(self, design: BinderDesign) -> Verdict:
        site = design.site
        ref = self.reference.geometry(site.metal)

        z = (site.bond_lengths() - ref.bond_length_mean) / ref.bond_length_std
        strain = float(np.sqrt(np.mean(z**2)))  # RMS bond-length deviation, in std units
        cn_gap = abs(site.coordination_number - ref.coordination_number)

        # Higher score = more plausible. Gaussian in geometric strain, penalised for
        # the wrong number of coordinating atoms.
        score = float(np.exp(-0.5 * strain**2) * np.exp(-cn_gap))
        ood = strain > self.ood_z
        trust = strain <= self.trust_z and cn_gap == 0 and not ood

        return Verdict(score=score, trust=trust, ood=ood, reason=self._reason(strain, cn_gap, ood))

    @staticmethod
    def _reason(strain: float, cn_gap: int, ood: bool) -> str:
        if ood:
            return f"off-manifold (bond strain {strain:.1f}σ) — defer"
        if cn_gap:
            return f"coordination number off by {cn_gap}"
        if strain > 2.0:
            return f"strained geometry ({strain:.1f}σ)"
        return f"plausible ({strain:.1f}σ, correct coordination)"
