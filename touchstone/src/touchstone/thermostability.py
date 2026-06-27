"""Whole-protein thermostability stage.

The MLIP-MD dynamics verifier checks whether the *coordination site* survives 300 K.
This stage asks the global question: does the whole protein stay folded — will it
hold up at the assay / operating temperature? A designed binder that melts at room
temperature is useless however good its site. This is the global counterpart to the
site-level MD survival check.

Pluggable predictor, mirroring the expression / co-fold stages: the heavy model
(a sequence Tm regressor like TemStaPro / DeepSTABp, or ThermoMPNN, or long MD) runs
on a GPU box (scripts/thermostability_score.py), so the verifier takes a
`predictor(design) -> ThermostabilitySignal | None`. `tm_provider` reads a precomputed
{sequence: Tm} map. Calibrate thresholds against measured Tm when available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .core import BinderDesign, Verdict, provider_from

_TM_SCALE = 10.0  # °C; fixes the score's sensitivity, independent of the trust cutoff


@dataclass
class ThermostabilitySignal:
    tm: float  # predicted melting temperature, °C


class ThermostabilityVerifier:
    """Trusts a design predicted to stay folded above `min_tm`; defers one predicted to
    be unfolded near room temperature (Tm below `ood_tm` — off the foldable manifold)
    or when no prediction is available."""

    def __init__(
        self,
        predictor: Callable[[BinderDesign], ThermostabilitySignal | None],
        *,
        min_tm: float = 50.0,
        ood_tm: float = 25.0,
    ):
        self.predictor = predictor
        self.min_tm = min_tm
        self.ood_tm = ood_tm

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            s = self.predictor(design)
        except Exception as e:  # model load / inference failure ⇒ can't judge
            return Verdict.defer(f"thermostability prediction failed: {type(e).__name__}")
        if s is None:
            return Verdict.defer("no thermostability prediction available")

        score = float(1.0 / (1.0 + np.exp(-(s.tm - self.min_tm) / _TM_SCALE)))  # 0.5 at min_tm
        reason = f"predicted Tm {s.tm:.0f} °C"
        if s.tm < self.ood_tm:  # unfolded near room temperature — off-manifold
            return Verdict.defer(reason, score=score)
        return Verdict(score, trust=s.tm >= self.min_tm, ood=False, reason=reason)


def tm_provider(predictions: dict[str, float]):
    """A predictor reading a precomputed {sequence: Tm °C} map (written by
    scripts/thermostability_score.py) — the in-library half of the stage."""
    return provider_from(predictions, transform=ThermostabilitySignal)
