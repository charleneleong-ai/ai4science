"""Angular coordination-geometry verifiers — the CheckMyMetal bases beyond bond length and
valence: nVECSUM (is the metal actually enclosed?) and polyhedron shape (is it the right
geometry, not just the right distances?). Both compute from the site alone (metal + donor
coordinates), so they run instantly with no extra input — closing the angular gap that
bond-length z-score + coordination number + bond-valence leave open.
"""

from __future__ import annotations

import numpy as np

from ..core import BinderDesign, Verdict


class CoordinationSymmetryVerifier:
    """Is the metal actually wrapped? The unit M→donor bond vectors of a complete, symmetric
    coordination sphere cancel — their mean (CheckMyMetal's nVECSUM, here the geometric
    unit-vector form) is ~0. A one-sided or incomplete site leaves a large residual: the metal
    is poking out, not enclosed. Trusts a balanced site; defers a lopsided one. This is
    invisible to bond-length / CN / valence checks — a site can have perfect distances and the
    right count yet leave the metal half-exposed."""

    def __init__(self, *, trust_v: float = 0.3, ood_v: float = 0.6):
        self.trust_v = trust_v  # |mean unit vector| within this ⇒ trusted (0 = balanced)
        self.ood_v = ood_v  # beyond this ⇒ one-sided ⇒ defer

    def verify(self, design: BinderDesign) -> Verdict:
        site = design.site
        if site.is_empty:
            return Verdict.defer("no coordinating atoms")
        v = site.ligand_xyz - site.metal_xyz
        vecsum = float(np.linalg.norm(np.mean(v / np.linalg.norm(v, axis=1, keepdims=True), axis=0)))
        score = float(np.exp(-0.5 * (vecsum / self.trust_v) ** 2))
        reason = f"coordination vector-sum {vecsum:.2f} (0 = enclosed, 1 = one-sided)"
        metrics = {"nvecsum": round(vecsum, 3), "cn": site.coordination_number}
        if vecsum > self.ood_v:
            return Verdict.defer(f"lopsided coordination (nVECSUM {vecsum:.2f})", score=score, metrics=metrics)
        return Verdict(score, trust=vecsum <= self.trust_v, ood=False, reason=reason, metrics=metrics)


# candidate ideal L–M–L angle multisets (degrees) per coordination number; each has
# C(cn,2) entries, matching `site.bond_angles()`. Multiple isomers per CN ⇒ best fit wins.
_IDEAL_ANGLES = {
    2: [[180.0]],
    3: [[120.0] * 3],
    4: [[109.47] * 6, [90.0] * 4 + [180.0] * 2],  # tetrahedral, square-planar
    5: [[90.0] * 6 + [120.0] * 3 + [180.0], [90.0] * 8 + [180.0] * 2],  # trig-bipyramidal, sq-pyramidal
    6: [[90.0] * 12 + [180.0] * 3],  # octahedral
}


class CoordinationGeometryVerifier:
    """Right distances, right shape? Scores the donor arrangement against the ideal
    coordination polyhedron for its CN — tetrahedral/square-planar at CN4,
    trigonal-bipyramidal/square-pyramidal at CN5, octahedral at CN6 — by RMS angle deviation
    (a gRMSD proxy, à la CheckMyMetal's geometry parameter). Catches a site with plausible
    bond lengths and the right count but a mangled, non-polyhedral arrangement."""

    def __init__(self, *, trust_deg: float = 20.0, ood_deg: float = 40.0):
        self.trust_deg = trust_deg  # RMS angle deviation within this ⇒ trusted
        self.ood_deg = ood_deg  # beyond this ⇒ off any ideal polyhedron ⇒ defer

    def verify(self, design: BinderDesign) -> Verdict:
        site = design.site
        if site.is_empty:
            return Verdict.defer("no coordinating atoms")
        ideals = _IDEAL_ANGLES.get(site.coordination_number)
        if not ideals:
            return Verdict.defer(f"no ideal polyhedron for CN={site.coordination_number}")
        observed = np.sort(site.bond_angles())
        rmsd = min(float(np.sqrt(np.mean((observed - np.sort(ideal)) ** 2))) for ideal in ideals)
        score = float(np.exp(-0.5 * (rmsd / self.trust_deg) ** 2))
        reason = f"polyhedron fit {rmsd:.1f}° RMS vs ideal CN{site.coordination_number}"
        metrics = {"angle_rmsd_deg": round(rmsd, 1), "cn": site.coordination_number}
        if rmsd > self.ood_deg:
            return Verdict.defer(f"distorted geometry ({rmsd:.0f}° off ideal)", score=score, metrics=metrics)
        return Verdict(score, trust=rmsd <= self.trust_deg, ood=False, reason=reason, metrics=metrics)
