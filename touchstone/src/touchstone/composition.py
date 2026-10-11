"""Model-free sequence sanity: reject degenerate amino-acid composition.

The one tier in the stack that runs no model at all, and it exists because the learned
ones are confidently wrong here: interface-confidence predictors score degenerate
sequences *highly*. In a 57-design conditional-binder campaign (tnf-binder, PR #81), 46
of the 53 designs that survived a memorisation check were low-complexity alanine/glycine
runs (59–76% alanine, one at 37% glycine), and one 59%-alanine sequence reached 0.556
binder-to-target ipTM — the top of a 40-design shortlist. A composition check is cheap,
model-free, and catches all of them immediately.

The expression tier is not a substitute: ESM-2 pseudo-perplexity is itself a learned
model, and a poly-alanine helix is not an improbable sequence to it.

Judges `design.sequence` only, so it stays orthogonal to the coordination site and to
every other tier. Opt-in rather than default-on: it judges the sequence, not the
coordination geometry the rest of the stack is calibrated against, so switching it on by
default would move every existing reward. Enable it with
`verify_structure(..., composition=True, sequence=...)` / `touchstone verify --composition
--sequence ...`, or hand it to `pipeline.cascade` / `pipeline.rank` directly.
"""

from __future__ import annotations

from collections import Counter

from .core import BinderDesign, Verdict

STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


class CompositionVerifier:
    """Trusts a sequence whose amino-acid usage resembles a real protein.

    Defaults are deliberately loose — they are meant to exclude degenerate output, not to
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

    @staticmethod
    def rejected(reason: str, metrics: dict | None = None) -> Verdict:
        """Judged and failed — weak, never defer. A degenerate composition is perfectly
        legible, it just isn't a protein's; so is a non-standard alphabet."""
        return Verdict(0.0, trust=False, ood=False, reason=reason, metrics=metrics or {})

    def verify(self, design: BinderDesign) -> Verdict:
        counts = Counter(design.sequence)
        if not counts:  # nothing to judge — distinct from judged-and-failed
            return Verdict.defer("empty sequence")
        if nonstandard := sorted(set(counts) - STANDARD_RESIDUES):
            return self.rejected(f"non-standard residues {nonstandard}")

        length = len(design.sequence)
        residue, count = counts.most_common(1)[0]
        top_fraction = count / length
        ala_gly = (counts["A"] + counts["G"]) / length
        metrics = {"top_residue": residue, "top_fraction": round(top_fraction, 3),
                   "distinct": len(counts), "ala_gly_fraction": round(ala_gly, 3)}

        if top_fraction > self.max_single_fraction:
            return self.rejected(f"low complexity: {residue} is {top_fraction:.0%}", metrics)
        if len(counts) < self.min_distinct:
            return self.rejected(f"only {len(counts)} distinct residues", metrics)
        if ala_gly > self.max_ala_gly_fraction:
            return self.rejected(f"Ala+Gly is {ala_gly:.0%} of the sequence", metrics)

        # mean of two already-normalised terms — headroom under the single-residue bound, and
        # alphabet coverage. Both land in [0, 1], so the score shares the scale the other tiers
        # report on: `reward_from_result` takes an fmean across them, and a half-scale tier
        # would quietly penalise every design that passes it.
        headroom = (self.max_single_fraction - top_fraction) / self.max_single_fraction
        score = (headroom + len(counts) / len(STANDARD_RESIDUES)) / 2
        return Verdict(round(score, 3), trust=True, ood=False, metrics=metrics,
                       reason=f"{len(counts)} distinct residues, max {residue} at {top_fraction:.0%}")
