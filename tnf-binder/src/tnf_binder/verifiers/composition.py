"""Reject degenerate sequences, independently of any model score.

This exists because interface-confidence models score poly-alanine designs confidently. One
such sequence (59% alanine) reached 0.556 binder-to-target ipTM and sat at the top of a
40-design shortlist. A composition check is cheap, model-free, and catches it immediately.
"""

from __future__ import annotations

from touchstone.core import Verdict

from tnf_binder.design import TnfDesign

STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


class CompositionVerifier:
    """Trusts a sequence whose amino-acid usage resembles a real protein.

    Defaults are deliberately loose -- they are meant to exclude degenerate output, not to
    impose a style on legitimate designs.
    """

    def __init__(
        self,
        *,
        max_single_fraction: float = 0.30,
        min_distinct: int = 12,
        max_ala_gly_fraction: float = 0.45,
    ) -> None:
        self.max_single_fraction = max_single_fraction
        self.min_distinct = min_distinct
        self.max_ala_gly_fraction = max_ala_gly_fraction

    def verify(self, design: TnfDesign) -> Verdict:
        counts = design.composition
        if not counts:
            return Verdict.defer("empty sequence")
        if design.has_nonstandard_residues:
            bad = sorted(set(design.sequence) - STANDARD_RESIDUES)
            return Verdict(0.0, trust=False, ood=False,
                           reason=f"non-standard residues {bad}")

        residue, count = counts.most_common(1)[0]
        top_fraction = count / design.length
        ala_gly = (counts["A"] + counts["G"]) / design.length
        metrics = {"top_residue": residue, "top_fraction": round(top_fraction, 3),
                   "distinct": len(counts), "ala_gly_fraction": round(ala_gly, 3)}

        if top_fraction >= self.max_single_fraction:
            return Verdict(0.0, trust=False, ood=False, metrics=metrics,
                           reason=f"low complexity: {residue} is {top_fraction:.0%}")
        if len(counts) < self.min_distinct:
            return Verdict(0.0, trust=False, ood=False, metrics=metrics,
                           reason=f"only {len(counts)} distinct residues")
        if ala_gly > self.max_ala_gly_fraction:
            return Verdict(0.0, trust=False, ood=False, metrics=metrics,
                           reason=f"Ala+Gly is {ala_gly:.0%} of the sequence")

        # mean of two already-normalised terms -- headroom under the single-residue bound, and
        # alphabet coverage. Clamping the sum to 1.0 before halving (the first version) made the
        # score exactly 0.500 for every real protein, since min_distinct=12 alone contributes
        # 0.6: a passing design could not be distinguished from any other passing design.
        headroom = (self.max_single_fraction - top_fraction) / self.max_single_fraction
        score = (headroom + len(counts) / len(STANDARD_RESIDUES)) / 2
        return Verdict(round(score, 3), trust=True, ood=False, metrics=metrics,
                       reason=f"{len(counts)} distinct residues, max {residue} at {top_fraction:.0%}")
