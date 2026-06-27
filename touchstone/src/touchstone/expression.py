"""Expression / sequence-plausibility stage.

The geometry, physics, and co-fold tiers all judge the *binding site*. This stage
asks the upstream wet-lab question: will the designed sequence actually express and
stay soluble? A binder that scores perfectly but won't come out of the cell is a
dead end. Combines two sequence-level signals:

  - ESM-2 pseudo-perplexity — how natural the sequence looks to a protein language
    model (lower = more expressible-looking).
  - a solubility score in [0, 1] (higher = more soluble).

The scorer is pluggable, mirroring the co-fold stage: the heavy ESM-2 forward +
solubility predictor run on a GPU box (scripts/expression_score.py), so the verifier
takes a `scorer(design) -> ExpressionSignals | None` — `None` when no precomputed
score is available (analogous to cofold's `provider`). `score_provider` reads a
precomputed {sequence: signals} map. Calibrate thresholds against wet-lab data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .core import BinderDesign, Verdict, reraise_if_bug


@dataclass
class ExpressionSignals:
    pseudo_perplexity: float  # ESM-2; lower = more natural
    solubility: float  # 0..1; higher = more soluble


class ExpressionVerifier:
    """Trusts a design whose sequence looks natural (low ESM-2 pseudo-perplexity)
    and soluble; defers a sequence so far off the natural manifold it can't be judged
    (pseudo-perplexity past `ood_perplexity`), or when no score is available."""

    def __init__(
        self,
        scorer: Callable[[BinderDesign], ExpressionSignals | None],
        *,
        max_perplexity: float = 12.0,
        ood_perplexity: float = 25.0,
        min_solubility: float = 0.4,
    ):
        self.scorer = scorer
        self.max_perplexity = max_perplexity
        self.ood_perplexity = ood_perplexity
        self.min_solubility = min_solubility

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            s = self.scorer(design)
        except Exception as e:  # model load / forward failure ⇒ can't judge
            reraise_if_bug(e)
            return Verdict.defer(f"expression scoring failed: {e}")
        if s is None:
            return Verdict.defer("no expression score available")

        # naturalness: 1 at/below the perplexity ceiling, decaying above it
        naturalness = float(np.clip(self.max_perplexity / max(s.pseudo_perplexity, 1e-6), 0.0, 1.0))
        score = naturalness * float(np.clip(s.solubility, 0.0, 1.0))
        reason = f"ESM-2 pseudo-perplexity {s.pseudo_perplexity:.1f}, solubility {s.solubility:.0%}"
        if s.pseudo_perplexity > self.ood_perplexity:  # off the natural-sequence manifold
            return Verdict.defer(reason, score=score)
        trust = s.pseudo_perplexity <= self.max_perplexity and s.solubility >= self.min_solubility
        return Verdict(score, trust=trust, ood=False, reason=reason)


def score_provider(scores: dict[str, ExpressionSignals]):
    """A scorer reading precomputed ExpressionSignals per design sequence (written by
    scripts/expression_score.py) — the in-library half of the expression stage."""

    def provide(design: BinderDesign) -> ExpressionSignals | None:
        return scores.get(design.sequence)

    return provide
