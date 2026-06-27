"""generator → verify → rank. The whole loop, in two functions."""

from __future__ import annotations

from dataclasses import replace

from .core import BinderDesign, Generator, Verdict, Verifier
from .geometry.ood import under_leachate, under_low_pH


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
