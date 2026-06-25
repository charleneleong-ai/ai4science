"""generator → verify → rank. The whole loop, in two functions."""

from __future__ import annotations

from .core import BinderDesign, Generator, Verdict, Verifier


def rank(designs: list[BinderDesign], verifier: Verifier) -> list[tuple[BinderDesign, Verdict]]:
    """Verify each design; return (design, verdict) sorted by score, best first."""
    scored = [(d, verifier.verify(d)) for d in designs]
    scored.sort(key=lambda dv: dv[1].score, reverse=True)
    return scored


def design_and_rank(
    generator: Generator, verifier: Verifier, target: str, n: int = 5
) -> list[tuple[BinderDesign, Verdict]]:
    return rank(generator.design(target, n), verifier)
