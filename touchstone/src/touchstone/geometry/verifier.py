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

        bonds = site.bond_lengths()
        if len(bonds) == 0:  # no coordinating atoms found — the worst case, not undefined
            return Verdict(0.0, trust=False, ood=True, reason="no coordinating atoms within cutoff — defer")

        z = (bonds - ref.bond_length_mean) / ref.bond_length_std
        strain = float(np.sqrt(np.mean(z**2)))  # RMS bond-length deviation, in std units
        cn_gap = abs(site.coordination_number - ref.coordination_number)  # vs modal, for the score
        cn_ok = ref.cn_range[0] <= site.coordination_number <= ref.cn_range[1]

        # Higher score = more plausible. Gaussian in geometric strain, penalised for
        # the distance from the modal coordination number.
        score = float(np.exp(-0.5 * strain**2) * np.exp(-cn_gap))
        ood = strain > self.ood_z
        trust = strain <= self.trust_z and cn_ok and not ood

        return Verdict(score=score, trust=trust, ood=ood, reason=self._reason(strain, cn_ok, ood))

    def _reason(self, strain: float, cn_ok: bool, ood: bool) -> str:
        if ood:
            return f"off-manifold (bond strain {strain:.1f}σ) — defer"
        if not cn_ok:
            return "coordination number outside observed range"
        if strain > self.trust_z:
            return f"strained geometry ({strain:.1f}σ)"
        return f"plausible ({strain:.1f}σ, coordination in range)"
