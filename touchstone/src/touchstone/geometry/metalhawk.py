"""MetalHawk geometry-distortion tier — an independent, licence-free geometry oracle.

MetalHawk (an artificial neural network trained on the CSD + MetalPDB —
github.com/vrettasm/MetalHawk) classifies a metal site's coordination number and
coordination geometry straight from the coordinates. It is the open stand-in for Mogul's
crystallographic geometry check: it corroborates the design's own CN and flags distorted
or ambiguous sites. Because it is a *learned classifier* (vs the geometry tier's empirical
z-score), agreement between the two is genuine defense-in-depth, not the same test twice.

Predictor-agnostic: the heavy ANN inference (running MetalHawk in its own venv) lives in
scripts/metalhawk_score.py and writes a prediction the verifier consumes via a `scorer`
callback — the same pluggable seam as the co-fold and expression tiers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from ..core import BinderDesign, Verdict, provider_from


@dataclass
class MetalHawkPrediction:
    coordination_number: int  # MetalHawk's predicted CN
    geometry: str  # predicted coordination-geometry class, e.g. "octahedral"
    confidence: float  # top-class softmax probability, 0..1


class MetalHawkVerifier:
    """Trusts a design whose coordination number MetalHawk independently predicts with high
    confidence; `weak` when the predicted CN disagrees (a distortion signal); `defer` when
    the prediction is too low-confidence to judge (off the ANN's training manifold) or no
    prediction is available. `scorer(design) -> MetalHawkPrediction | None`."""

    def __init__(self, scorer: Callable[[BinderDesign], MetalHawkPrediction | None],
                 *, min_confidence: float = 0.5, ood_confidence: float = 0.3):
        self.scorer = scorer
        self.min_confidence = min_confidence
        self.ood_confidence = ood_confidence

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            p = self.scorer(design)
        except Exception as e:  # inference / parse failure ⇒ can't judge
            return Verdict.defer(f"MetalHawk scoring failed: {type(e).__name__}")
        if p is None:
            return Verdict.defer("no MetalHawk prediction available")

        cn_design = design.site.coordination_number
        cn_match = p.coordination_number == cn_design
        score = float(p.confidence * (1.0 if cn_match else 0.5))  # confident agreement scores highest
        metrics = {"cn_pred": p.coordination_number, "cn_design": cn_design,
                   "geometry": p.geometry, "confidence": round(p.confidence, 3)}
        reason = (f"MetalHawk {p.geometry} CN{p.coordination_number} vs design CN{cn_design}, "
                  f"confidence {p.confidence:.0%}")
        if p.confidence < self.ood_confidence:  # ANN unsure ⇒ off the training manifold, not judgeable
            return Verdict.defer(reason, score=score, metrics=metrics)
        trust = cn_match and p.confidence >= self.min_confidence
        return Verdict(score, trust=trust, ood=False, reason=reason, metrics=metrics)


def score_provider(scores: dict[str, MetalHawkPrediction]):
    """A scorer reading precomputed MetalHawk predictions keyed by structure path (written by
    scripts/metalhawk_score.py) — the in-library half of the MetalHawk tier."""
    return provider_from(scores, key="source")
