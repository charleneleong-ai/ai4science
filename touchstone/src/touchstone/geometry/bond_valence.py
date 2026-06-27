"""A second geometry verifier, by a different principle than the z-score oracle.

The bond-valence model (Brown–Altermatt; Brese & O'Keeffe) says each metal–donor
bond carries a valence vᵢ = exp((R₀ − dᵢ)/b), and that the sum over a site's bonds
should recover the metal's formal oxidation state. A bond-valence sum (BVS) far
from the formal charge flags an over- or under-bonded — geometrically implausible —
site, judged independently of how the distances compare to a reference distribution.
Pairs with GeometryVerifier the way two independent co-folders or two MLIPs do.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from ..core import BinderDesign, Verdict, oxidation_state

_PARAMS = Path(__file__).parent.parent / "data" / "bond_valence_params.json"


def _load_params() -> dict:
    return json.loads(_PARAMS.read_text())


class BondValenceVerifier:
    """Trusts a site whose bond-valence sum recovers the metal's formal charge.

    `trust_tol` / `ood_tol` are valence-unit deviations of the BVS from the formal
    oxidation state (e.g. 2 for Ni²⁺).
    """

    def __init__(self, params: dict | None = None, trust_tol: float = 0.4, ood_tol: float = 0.8):
        self.params = params or _load_params()
        self.trust_tol = trust_tol
        self.ood_tol = ood_tol

    def verify(self, design: BinderDesign) -> Verdict:
        site = design.site
        r0 = self.params["metals"][site.metal]  # KeyError on unknown metal — like GeometryVerifier
        b = self.params["b"]
        valence = oxidation_state(site.metal)

        if site.is_empty:
            return Verdict.defer("no coordinating atoms")
        # Don't silently undercount: an unparameterized donor would drop out of the sum
        # and wrongly deflate the BVS — defer rather than judge a partial valence.
        unparameterized = sorted(set(site.ligand_elems) - set(r0))
        if unparameterized:
            return Verdict.defer(f"donor element(s) {','.join(unparameterized)} not parameterized")

        bvs = sum(math.exp((r0[e] - d) / b) for e, d in zip(site.ligand_elems, site.bond_lengths()))
        disc = abs(bvs - valence)
        score = math.exp(-0.5 * (disc / self.trust_tol) ** 2)
        reason = f"BVS {bvs:.2f} vs formal {valence} (Δ{disc:.2f})"
        if disc > self.ood_tol:
            return Verdict.defer(reason, score=score)
        return Verdict(score, trust=disc <= self.trust_tol, ood=False, reason=reason)
