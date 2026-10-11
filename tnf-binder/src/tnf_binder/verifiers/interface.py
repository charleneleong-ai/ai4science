"""Judge the designed interface from buried surface area, not contact counts.

Counting atom or residue pairs was tried first and was actively misleading: the count scales
with interpenetration, so clashing poly-alanine designs out-scored real ones. Buried surface
area (unbound SASA minus bound SASA) is a direct physical quantity and degrades gracefully.

Separately: pre-fold geometry does NOT predict interface confidence. On folded designs the
cleanest geometry ranked 5th of 6 while flagged junk ranked 1st, so this verifier judges a
structure, and never selects which structures are worth folding.
"""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import dataclass

from touchstone.core import Verdict

from tnf_binder.design import TnfDesign

MIN_INTERFACE_AREA = 350.0
MIN_EPITOPE_FOCUS = 0.25


@dataclass(frozen=True)
class InterfaceMetrics:
    """Buried-area view of one binder/target interface."""

    area: float                 # target-side buried area, A^2
    epitope_focus: float        # share of that area on the intended epitope
    divergent_fraction: float   # share on human/mouse-divergent positions
    clashes: int

    @classmethod
    def from_areas(
        cls,
        areas: dict[int, float],
        *,
        epitope: Collection[int],
        divergent: Collection[int],
        clashes: int,
    ) -> "InterfaceMetrics":
        """Summarise per-position buried area against the intended epitope.

        Shares are of buried area, not residue counts: one deeply-buried residue matters more
        than three glancing ones, and counting treats them alike.
        """
        total = sum(areas.values())

        def share(positions: Collection[int]) -> float:
            if not total:
                return 0.0
            return sum(area for p, area in areas.items() if p in positions) / total

        return cls(area=total, epitope_focus=share(epitope),
                   divergent_fraction=share(divergent), clashes=clashes)

    @property
    def cross_reactivity(self) -> float:
        """Objective 2 proxy: the share of interface NOT on divergent residues."""
        return 1.0 - self.divergent_fraction


class InterfaceVerifier:
    """Trusts a real interface, on the intended epitope, without steric clashes.

    A real minibinder interface buries roughly 600-1000 A^2 per side.
    """

    def __init__(
        self,
        *,
        min_area: float = MIN_INTERFACE_AREA,
        min_epitope_focus: float = MIN_EPITOPE_FOCUS,
        max_clashes: int = 0,
    ) -> None:
        self.min_area = min_area
        self.min_epitope_focus = min_epitope_focus
        self.max_clashes = max_clashes

    def verify(self, design: TnfDesign) -> Verdict:
        metrics: InterfaceMetrics | None = design.metrics.get("interface")
        if metrics is None:
            return Verdict.defer("no interface measurement available")

        detail = {
            "area": round(metrics.area, 1),
            "epitope_focus": round(metrics.epitope_focus, 3),
            "cross_reactivity": round(metrics.cross_reactivity, 3),
            "clashes": metrics.clashes,
        }
        if metrics.clashes > self.max_clashes:
            return Verdict(0.0, trust=False, ood=False, metrics=detail,
                           reason=f"{metrics.clashes} steric clashes")
        if metrics.area < self.min_area:
            return Verdict(0.0, trust=False, ood=False, metrics=detail,
                           reason=f"interface only {metrics.area:.0f} A^2")
        if metrics.epitope_focus < self.min_epitope_focus:
            return Verdict(0.0, trust=False, ood=False, metrics=detail,
                           reason=f"off epitope: {metrics.epitope_focus:.0%} of the interface")

        score = min(1.0, metrics.area / 1000.0) * metrics.epitope_focus
        return Verdict(round(score, 3), trust=True, ood=False, metrics=detail,
                       reason=f"{metrics.area:.0f} A^2 buried, {metrics.epitope_focus:.0%} on "
                              f"epitope, {metrics.cross_reactivity:.0%} conserved")
