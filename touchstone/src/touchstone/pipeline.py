"""generator → verify → rank, and a cost-ordered cascade over the tiers."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

from .core import BinderDesign, Generator, Verdict, Verifier
from .geometry.ood import under_leachate, under_low_pH


@dataclass
class CascadeResult:
    """How one design fared through a cost-ordered cascade of verifier tiers."""

    design: BinderDesign
    verdicts: list[tuple[str, Verdict]]  # (tier label, verdict), in the order run
    survived: bool  # advanced past every tier

    @property
    def dropped_at(self) -> str | None:
        """Tier that rejected the design, or None if it survived all of them."""
        return None if self.survived else self.verdicts[-1][0]

    @property
    def final(self) -> Verdict:
        return self.verdicts[-1][1]


def _advances(verdict: Verdict) -> bool:
    """Default gate: a design advances to the next tier unless flagged off-manifold."""
    return not verdict.ood


def cascade(
    designs: list[BinderDesign],
    tiers: list[tuple[str, Verifier]],
    *,
    advances: Callable[[Verdict], bool] = _advances,
) -> list[CascadeResult]:
    """Run verifier `tiers` (label, verifier) cheap→expensive; each design advances only
    while `advances(verdict)` holds, stopping at the first tier that rejects it — so the
    expensive tiers (MLIP / co-fold / MD) only ever run on designs the cheap CSD/geometry
    gates let through. Returns a CascadeResult per design.

    `advances` defaults to "not off-manifold" (drop only clear rejects); pass
    `lambda v: v.trust` for aggressive triage that only fast-tracks trusted designs.
    """
    results = []
    for design in designs:
        verdicts: list[tuple[str, Verdict]] = []
        survived = True
        for label, verifier in tiers:
            verdict = verifier.verify(design)
            verdicts.append((label, verdict))
            if not advances(verdict):
                survived = False
                break
        results.append(CascadeResult(design, verdicts, survived))
    return results


def rank(designs: list[BinderDesign], verifier: Verifier) -> list[tuple[BinderDesign, Verdict]]:
    """Verify each design; return (design, verdict) sorted by score, best first."""
    scored = [(d, verifier.verify(d)) for d in designs]
    scored.sort(key=lambda dv: dv[1].score, reverse=True)
    return scored


def design_and_rank(
    generator: Generator, verifier: Verifier, target: str, n: int = 5
) -> list[tuple[BinderDesign, Verdict]]:
    return rank(generator.design(target, n), verifier)


def selectivity_profile(
    design: BinderDesign, verifier: Verifier, metals: list[str]
) -> dict[str, Verdict]:
    """Re-score a design's site as each competing metal → {metal: Verdict}.

    A genuinely selective binder trusts for its target metal and not the competitors.
    Caveat the verdict reveals: geometry alone rarely discriminates divalent metals —
    their coordination distributions overlap, so this usually returns `trust` for all of
    them. Real selectivity needs donor-identity (HSAB) or physics (differential binding),
    not coordination geometry.
    """
    return {m: verifier.verify(replace(design, site=replace(design.site, metal=m))) for m in metals}


def stress_profile(
    design: BinderDesign, verifier: Verifier, *, bond_stretch: float = 0.5, n_protonate: int = 1
) -> dict[str, Verdict]:
    """Re-verify a design under extreme-condition perturbations → {condition: Verdict}.

    A robust binder holds its verdict across the operating envelope; a fragile one
    degrades to weak/defer under stress. Conditions: `neutral` (as-is), `leachate`
    (bonds stretched — hot/acidic/saline), `low_pH` (the most labile donors protonated
    off). The robustness map, not a single verdict — does it work in the real process?
    """
    return {
        "neutral": verifier.verify(design),
        "leachate": verifier.verify(under_leachate(design, bond_stretch)),
        "low_pH": verifier.verify(under_low_pH(design, n_protonate)),
    }
