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

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from .core import BinderDesign, Verdict, provider_from


@dataclass
class ExpressionSignals:
    pseudo_perplexity: float  # ESM-2; lower = more natural
    solubility: float  # 0..1; higher = more soluble


class ExpressionVerifier:
    """Rejects a sequence far enough off the natural manifold that it usually fails to express,
    and declines to rank anything below that; defers past `ood_perplexity` or with no score.

    **Perplexity is a junk filter, not a ranker, and this class used to treat it as one.** The
    score was `max_perplexity / perplexity`, a monotonic decrease, so a design at perplexity 5
    outscored one at 10. Measured against 681 designs with an experimental expression verdict
    -- every non-expressing design in the Proteinbase export, length-matched to a success --
    that gradient is not there and locally runs the other way:

        perplexity   0-5    5-7    7-9   9-12  12-16    16+
        expressed  45.8%  41.8%  63.0%  70.6%  44.2%   9.3%

    Expression *peaks* at 9-12, the region the old score penalised most. The pattern is not
    length in disguise: mean length is flat across those bins (113-158 residues) and the
    non-monotonicity survives within a single length band. What does hold is the tail -- above
    perplexity 16, expression is 9.3% against 52.7%, a plateau whose separation clears
    Bonferroni over 33 scanned cuts from 12 through 20.

    So the score is now flat below the cut rather than graded, because a gradient the data does
    not show is worse than no gradient: it ranks designs confidently in the wrong order.
    `max_perplexity` moves 12.0 -> 16.0 to sit at the measured plateau. Both are defaults from
    one export pooling many targets and expression systems; re-run `touchstone calibrate`
    against your own outcomes before relying on either.

    `min_solubility` is **not** calibrated and cannot be from this data: the solubility signal
    is the placeholder GRAVY heuristic, and expression is flat in GRAVY across the whole export.
    """

    def __init__(
        self,
        scorer: Callable[[BinderDesign], ExpressionSignals | None],
        *,
        max_perplexity: float = 16.0,
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
            return Verdict.defer(f"expression scoring failed: {type(e).__name__}")
        if s is None:
            return Verdict.defer("no expression score available")

        # Flat below the cut, falling away above it. Below the cut the measurement shows no
        # usable gradient, so scoring one differently from another would be invention.
        over = s.pseudo_perplexity > self.max_perplexity
        naturalness = (float(np.clip(self.max_perplexity / s.pseudo_perplexity, 0.0, 1.0))
                       if over else 1.0)
        score = naturalness * float(np.clip(s.solubility, 0.0, 1.0))
        reason = f"ESM-2 pseudo-perplexity {s.pseudo_perplexity:.1f}, solubility {s.solubility:.0%}"
        if s.pseudo_perplexity > self.ood_perplexity:  # off the natural-sequence manifold
            return Verdict.defer(reason, score=score)
        trust = s.pseudo_perplexity <= self.max_perplexity and s.solubility >= self.min_solubility
        return Verdict(score, trust=trust, ood=False, reason=reason)


def score_provider(scores: dict[str, ExpressionSignals]):
    """A scorer reading precomputed ExpressionSignals per design sequence (written by
    scripts/expression_score.py) — the in-library half of the expression stage."""
    return provider_from(scores)


def load_signals(path: str | Path) -> dict[str, ExpressionSignals]:
    """Load an expression-score JSON ({sequence: {pseudo_perplexity, solubility}}) into
    {sequence: ExpressionSignals} for score_provider — the loader that closes the file→verifier loop."""
    raw = json.loads(Path(path).read_text())
    return {seq: ExpressionSignals(**v) for seq, v in raw.items()}
