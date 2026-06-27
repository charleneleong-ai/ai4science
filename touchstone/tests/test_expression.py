import numpy as np
import pytest

from touchstone import BinderDesign, ExpressionSignals, ExpressionVerifier
from touchstone.core import CoordinationSite


def _design(seq: str = "MKLVINGKTLKGEITV") -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
    return BinderDesign(seq, site, generator="test", generator_confidence=0.5)


def _scorer(ppl: float, sol: float):
    def score(_design: BinderDesign) -> ExpressionSignals:
        return ExpressionSignals(pseudo_perplexity=ppl, solubility=sol)

    return score


class TestExpressionVerifier:
    def test_natural_and_soluble_trusts(self):
        v = ExpressionVerifier(_scorer(ppl=6.0, sol=0.8)).verify(_design())
        assert v.trust and v.score > 0.5

    @pytest.mark.parametrize("ppl, sol", [(18.0, 0.8), (6.0, 0.1)], ids=["unnatural", "insoluble"])
    def test_failing_either_signal_is_weak(self, ppl, sol):
        # judgeable but not trusted — neither off-manifold nor a failure
        v = ExpressionVerifier(_scorer(ppl=ppl, sol=sol)).verify(_design())
        assert not v.trust and not v.ood

    def test_off_manifold_sequence_defers(self):
        v = ExpressionVerifier(_scorer(ppl=40.0, sol=0.8)).verify(_design())  # past ood_perplexity
        assert v.ood and not v.trust

    def test_score_rewards_both_signals(self):
        good = ExpressionVerifier(_scorer(6.0, 0.9)).verify(_design()).score
        meh = ExpressionVerifier(_scorer(10.0, 0.5)).verify(_design()).score
        assert good > meh

    def test_no_score_available_defers(self):
        v = ExpressionVerifier(lambda _d: None).verify(_design())
        assert v.ood and "no expression score" in v.reason

    def test_scorer_failure_defers(self):
        def boom(_design):
            raise RuntimeError("model died")

        v = ExpressionVerifier(boom).verify(_design())
        assert v.ood and not v.trust and "failed" in v.reason
