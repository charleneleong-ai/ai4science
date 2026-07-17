"""A second geometry verifier, by a different principle than the z-score oracle.

The bond-valence model (Brown–Altermatt; Brese & O'Keeffe) says each metal–donor
bond carries a valence vᵢ = exp((R₀ − dᵢ)/b), and that the sum over a site's bonds
should recover the metal's formal oxidation state. A bond-valence sum (BVS) far
from the formal charge flags an over- or under-bonded — geometrically implausible —
site, judged independently of how the distances compare to a reference distribution.
Pairs with GeometryVerifier the way two independent co-folders or two MLIPs do.

The R₀ and tolerances are **metalloprotein-domain** values (scripts/build_bond_valence_params.py),
not the small-molecule Brese & O'Keeffe originals: those under-count protein sites (median real Ni²⁺
BVS ≈ 1.03 vs formal 2), the same wrong-domain error as the geometry prior. Recalibrated, BVS is a
**weak** filter in this domain — protein coordination is heterogeneous, so the sum's scatter is wide
(wider for Ni/Co than Cu) — hence per-metal tolerances read from the reference data. Its role is
catching egregious over/under-bonding, not fine discrimination; the ensemble discriminates.
See docs/experiments/2026-07-17-bond-valence-wrong-domain.md.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from ..core import BinderDesign, Verdict, oxidation_state

PARAMS = Path(__file__).parent.parent / "data" / "bond_valence_params.json"
TOLERANCE_KEYS = ("trust_tol", "ood_tol", "domain")  # per-metal metadata, not donor R₀ values


def load_params() -> dict:
    return json.loads(PARAMS.read_text())


class BondValenceVerifier:
    """Trusts a site whose bond-valence sum recovers the metal's formal charge.

    Tolerances are the valence-unit deviation of the BVS from the formal oxidation state at which a
    site stops trusting (`trust_tol`) or is deferred as off-manifold (`ood_tol`). They default to the
    per-metal values in the params (the protein-domain BVS scatter); the constructor args are only a
    fallback for a metal whose params carry none, or for an injected test table."""

    def __init__(self, params: dict | None = None, trust_tol: float = 0.4, ood_tol: float = 0.8):
        self.params = params or load_params()
        self.trust_tol = trust_tol
        self.ood_tol = ood_tol

    def verify(self, design: BinderDesign) -> Verdict:
        site = design.site
        entry = self.params["metals"][site.metal]  # KeyError on unknown metal — like GeometryVerifier
        r0 = {e: v for e, v in entry.items() if e not in TOLERANCE_KEYS}
        trust_tol = entry.get("trust_tol", self.trust_tol)  # data owns the tolerance; ctor is fallback
        ood_tol = entry.get("ood_tol", self.ood_tol)
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
        score = math.exp(-0.5 * (disc / trust_tol) ** 2)
        reason = f"BVS {bvs:.2f} vs formal {valence} (Δ{disc:.2f})"
        metrics = {"bvs": round(bvs, 2), "formal_valence": valence, "delta": round(disc, 2)}
        if disc > ood_tol:
            return Verdict.defer(reason, score=score, metrics=metrics)
        return Verdict(score, trust=disc <= trust_tol, ood=False, reason=reason, metrics=metrics)
