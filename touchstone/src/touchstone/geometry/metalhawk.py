"""MetalHawk geometry-distortion tier — EXPERIMENTAL, demoted (open ANN geometry oracle).

**Status: experimental, not for de-novo designs.** End-to-end evaluation on 96 BoltzGen
designs (docs/experiments/2026-07-07-metalhawk-ood-designed-sites.md) showed MetalHawk is
*confidently out-of-distribution* on designed sites — it calls 95/96 clean CN4 sites
octahedral (CN6) at ~100% confidence. Wired into the consensus it collapsed 36 trusts → 0.
The analytic `coord_geometry` (polyhedron-RMSD) tier is the geometry oracle for designs; this
tier is a licence-free option only for in-distribution (natural-like) sites. Unlike Mogul —
an empirical CSD *lookup* that abstains ("insufficient data") when off-manifold — a learned
classifier extrapolates confidently, which is exactly the failure mode here.

MetalHawk (an artificial neural network trained on the CSD + MetalPDB —
github.com/vrettasm/MetalHawk) classifies a metal site's coordination geometry (LIN/TRI/
TET/SPL/SQP/TBP/OCT, and thus CN) straight from the coordinates. It is the open, learned
stand-in for Mogul's crystallographic geometry check — a second opinion on the design's CN.

Because it is a *learned* classifier, it goes out-of-distribution: on de-novo designed
sites (off its CSD/MetalPDB manifold) it can be *confidently wrong* — e.g. calling a clean
4-coordinate site octahedral (CN6) at ~100% confidence. So this tier trusts MetalHawk only
when it agrees with the structure's own distance-derived CN, and — crucially — treats a
*confident yet gross* CN disagreement (Δcn ≥ `ood_cn_gap`) as an off-manifold signal and
`defer`s, rather than emitting a spurious `weak`. The deterministic `coord_geometry`
polyhedron-RMSD tier carries the real geometry signal; MetalHawk is an optional cross-check
that abstains when out of its depth. A small Δcn (a genuine mild distortion) is still `weak`.

Predictor-agnostic: the heavy ANN inference (running MetalHawk in its own venv) lives in
scripts/metalhawk_score.py and writes a JSON the verifier consumes via `load_predictions`
+ `score_provider` — the same pluggable file-handoff seam as the co-fold / expression tiers.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from ..core import BinderDesign, Verdict, provider_from


@dataclass
class MetalHawkPrediction:
    coordination_number: int  # CN implied by MetalHawk's geometry class (OCT→6, TET→4, …)
    geometry: str  # predicted coordination-geometry class, e.g. "octahedral"
    confidence: float  # 1 − entropy/ln(7), MetalHawk's inverse-entropy certainty, 0..1


class MetalHawkVerifier:
    """`trust` when MetalHawk's CN matches the structure's own distance-derived CN with high
    confidence; `weak` on a small confident mismatch (a mild distortion signal); `defer` when
    MetalHawk is unsure (< `ood_confidence`) **or** confidently contradicts the physical shell
    by Δcn ≥ `ood_cn_gap` (confidently off its training manifold — can't be trusted as signal);
    `defer` too when no prediction is available. `scorer(design) -> MetalHawkPrediction | None`."""

    def __init__(self, scorer: Callable[[BinderDesign], MetalHawkPrediction | None],
                 *, min_confidence: float = 0.5, ood_confidence: float = 0.3, ood_cn_gap: int = 2):
        self.scorer = scorer
        self.min_confidence = min_confidence
        self.ood_confidence = ood_confidence
        self.ood_cn_gap = ood_cn_gap

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            p = self.scorer(design)
        except Exception as e:  # inference / parse failure ⇒ can't judge
            return Verdict.defer(f"MetalHawk scoring failed: {type(e).__name__}")
        if p is None:
            return Verdict.defer("no MetalHawk prediction available")

        cn_design = design.site.coordination_number  # the structure's own distance-derived CN
        cn_gap = abs(p.coordination_number - cn_design)
        score = float(p.confidence * (1.0 if cn_gap == 0 else 0.5))  # confident agreement scores highest
        metrics = {"cn_pred": p.coordination_number, "cn_design": cn_design, "cn_gap": cn_gap,
                   "geometry": p.geometry, "confidence": round(p.confidence, 3)}
        reason = (f"MetalHawk {p.geometry} CN{p.coordination_number} vs design CN{cn_design}, "
                  f"confidence {p.confidence:.0%}")
        if p.confidence < self.ood_confidence:  # the ANN itself is unsure ⇒ off-manifold, not judgeable
            return Verdict.defer(f"{reason} — low confidence", score=score, metrics=metrics)
        if cn_gap >= self.ood_cn_gap:  # confident yet grossly wrong on the physical CN ⇒ confidently off-manifold
            return Verdict.defer(f"{reason}; off-manifold (Δcn {cn_gap} vs physical shell)",
                                 score=score, metrics=metrics)
        trust = cn_gap == 0 and p.confidence >= self.min_confidence
        return Verdict(score, trust=trust, ood=False, reason=reason, metrics=metrics)


def load_predictions(path: str | Path) -> dict[str, MetalHawkPrediction]:
    """Load scripts/metalhawk_score.py's JSON into {structure_path: MetalHawkPrediction},
    ready for score_provider — the loader that closes the file-handoff loop."""
    raw = json.loads(Path(path).read_text())
    return {k: MetalHawkPrediction(v["coordination_number"], v["geometry"], v["confidence"])
            for k, v in raw.items()}


def score_provider(scores: dict[str, MetalHawkPrediction]):
    """A scorer reading precomputed MetalHawk predictions keyed by structure path (written by
    scripts/metalhawk_score.py, loaded via `load_predictions`) — the in-library half of the tier."""
    return provider_from(scores, key="source")
