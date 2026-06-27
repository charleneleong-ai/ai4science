import importlib.util
from pathlib import Path

import numpy as np
import pytest

from touchstone import BinderDesign, ThermostabilitySignal, ThermostabilityVerifier
from touchstone.core import CoordinationSite

# the TemStaPro scorer runs in its own conda env so it can't import touchstone; load it by
# path to unit-test the pure threshold→Tm mapping (the contract between the model and the verifier)
_spec = importlib.util.spec_from_file_location(
    "thermostability_score", Path(__file__).parent.parent / "scripts" / "thermostability_score.py"
)
tsp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tsp)


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


class TestTemStaProScoring:
    """The TemStaPro classifier→Tm summarisation feeding tm_provider."""

    def _passed(self, *up_to: int) -> dict[int, bool]:
        return {t: t in up_to for t in tsp.THRESHOLDS}

    @pytest.mark.parametrize(
        "up_to, expected",
        [
            ((), 37.0),  # not stable above even 40 °C ⇒ mesophilic/unstable
            ((40,), 40.0),  # stable above 40 only
            ((40, 45, 50), 50.0),  # highest passed threshold drives the estimate
            (tsp.THRESHOLDS, 70.0),  # stable through 65 ⇒ thermophilic, saturates
        ],
    )
    def test_tm_is_highest_passed_threshold(self, up_to, expected):
        assert tsp._tm_from_thresholds(self._passed(*up_to)) == expected

    def test_binary_columns_prefers_exact_over_raw(self):
        # mean-output carries a binary + a 'raw' probability column per threshold; pick the binary
        fields = ["sequence"] + [f"{t}" for t in tsp.THRESHOLDS] + [f"{t}_raw" for t in tsp.THRESHOLDS] + ["clash"]
        cols = tsp._binary_columns(fields)
        assert cols == {t: str(t) for t in tsp.THRESHOLDS}
