import numpy as np
import pytest

from touchstone import BinderDesign, ThermostabilitySignal, ThermostabilityVerifier
from touchstone.core import CoordinationSite


def _design(seq: str = "MKLVINGKTLKGEITV") -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
    return BinderDesign(seq, site, generator="test", generator_confidence=0.5)


def _pred(tm: float):
    def predict(_design: BinderDesign) -> ThermostabilitySignal:
        return ThermostabilitySignal(tm=tm)

    return predict


class TestThermostabilityVerifier:
    def test_stable_design_trusts(self):
        v = ThermostabilityVerifier(_pred(70.0)).verify(_design())
        assert v.trust and v.score > 0.5

    def test_marginal_design_is_weak(self):
        # 25 °C (ood) < 40 < 50 (min): judgeable but not trusted
        v = ThermostabilityVerifier(_pred(40.0)).verify(_design())
        assert not v.trust and not v.ood

    def test_room_temp_unfolded_defers(self):
        v = ThermostabilityVerifier(_pred(15.0)).verify(_design())  # below ood_tm ⇒ off-manifold
        assert v.ood and not v.trust

    def test_higher_tm_scores_higher(self):
        hi = ThermostabilityVerifier(_pred(80.0)).verify(_design()).score
        lo = ThermostabilityVerifier(_pred(55.0)).verify(_design()).score
        assert hi > lo

    def test_no_prediction_defers(self):
        v = ThermostabilityVerifier(lambda _d: None).verify(_design())
        assert v.ood and "no thermostability" in v.reason

    def test_prediction_failure_defers(self):
        def boom(_design):
            raise RuntimeError("model died")

        v = ThermostabilityVerifier(boom).verify(_design())
        assert v.ood and not v.trust and "failed" in v.reason
