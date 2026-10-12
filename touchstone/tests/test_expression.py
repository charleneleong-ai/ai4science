import json

import numpy as np
import pytest

from touchstone import BinderDesign, ExpressionSignals, ExpressionVerifier, load_signals
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

    def test_solubility_still_grades_the_score(self):
        """It used to be perplexity doing this work; now only solubility varies below the cut."""
        good = ExpressionVerifier(_scorer(6.0, 0.9)).verify(_design()).score
        meh = ExpressionVerifier(_scorer(6.0, 0.5)).verify(_design()).score
        assert good > meh

    @pytest.mark.parametrize("low,high", [(5.0, 10.0), (6.0, 15.9), (1.0, 16.0)])
    def test_perplexity_does_not_grade_below_the_cut(self, low, high):
        """Measured over 681 designs, expression is 45.8% at perplexity 0-5 and 70.6% at 9-12,
        so the old `max_perplexity / perplexity` gradient ranked designs confidently backwards.
        Flat is what the data supports."""
        scores = [ExpressionVerifier(_scorer(p, 0.8)).verify(_design()).score for p in (low, high)]
        assert scores[0] == pytest.approx(scores[1])

    def test_the_cut_sits_at_the_measured_plateau(self):
        """Above perplexity 16 expression is 9.3% against 52.7%; 14 is inside the flat region
        and used to be rejected by the old default of 12."""
        assert ExpressionVerifier(_scorer(14.0, 0.8)).verify(_design()).trust
        assert not ExpressionVerifier(_scorer(17.0, 0.8)).verify(_design()).trust

    def test_the_score_still_falls_away_above_the_cut(self):
        over = ExpressionVerifier(_scorer(20.0, 0.8)).verify(_design()).score
        under = ExpressionVerifier(_scorer(10.0, 0.8)).verify(_design()).score
        assert over < under

    def test_no_score_available_defers(self):
        v = ExpressionVerifier(lambda _d: None).verify(_design())
        assert v.ood and "no expression score" in v.reason

    def test_scorer_failure_defers(self):
        def boom(_design):
            raise RuntimeError("model died")

        v = ExpressionVerifier(boom).verify(_design())
        assert v.ood and not v.trust and "failed" in v.reason


def test_load_signals_round_trips(tmp_path):
    # the loader that closes the expression-score JSON → verifier loop
    p = tmp_path / "expr.json"
    p.write_text(json.dumps({"SEQ": {"pseudo_perplexity": 7.0, "solubility": 0.8}}))
    assert load_signals(p) == {"SEQ": ExpressionSignals(7.0, 0.8)}
