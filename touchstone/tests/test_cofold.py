import numpy as np
import pytest

from touchstone import BinderDesign, CofoldCrossCheck, cofold_agreement
from touchstone.core import CoordinationSite

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def _site(elems=("N", "N", "O", "O", "N", "O"), bond: float = 2.1) -> CoordinationSite:
    n = len(elems)
    return CoordinationSite("Ni2+", np.zeros(3), _OCT[:n] * bond, tuple(elems))


def _design(site: CoordinationSite) -> BinderDesign:
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5, source="d.cif")


def _fixed(site):
    """A provider returning `site` for any design — a stand-in for a real predictor."""

    def provide(_design):
        return site

    return provide


class TestCofoldAgreement:
    def test_matching_prediction_trusts(self):
        v = cofold_agreement(_site(), _site())  # identical site recovered
        assert v.trust and v.score > 0.9

    def test_cn_mismatch_is_weak_not_deferred(self):
        # an independent predictor placing a very different CN ⇒ not corroborated,
        # but still a judgeable verdict (weak), not defer
        v = cofold_agreement(_site(), _site(elems=("N", "O")))
        assert not v.trust and not v.ood

    def test_wrong_donor_identity_lowers_score(self):
        # design wants all-N coordination; predictor places all-O ⇒ no donor overlap
        v = cofold_agreement(_site(elems=("N",) * 6), _site(elems=("O",) * 6))
        assert not v.trust and v.score < 0.5

    def test_stretched_prediction_not_trusted(self):
        v = cofold_agreement(_site(bond=2.1), _site(bond=2.9))  # Δbond 0.8 > tol
        assert not v.trust

    def test_empty_prediction_defers(self):
        empty = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
        v = cofold_agreement(_site(), empty)
        assert v.ood and not v.trust

    def test_empty_reference_defers_with_finite_score(self):
        # design's own site empty ⇒ defer, not a NaN score (which would corrupt rank())
        empty = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
        v = cofold_agreement(empty, _site())
        assert v.ood and not v.trust and np.isfinite(v.score)


class TestCofoldCrossCheck:
    def test_agreeing_provider_trusts(self):
        v = CofoldCrossCheck(provider=_fixed(_site())).verify(_design(_site()))
        assert v.trust

    def test_missing_prediction_defers(self):
        v = CofoldCrossCheck(provider=_fixed(None)).verify(_design(_site()))
        assert v.ood and "no co-fold" in v.reason

    def test_provider_error_defers(self):
        def boom(_):
            raise RuntimeError("inference died")

        v = CofoldCrossCheck(provider=boom).verify(_design(_site()))
        assert v.ood and "inference died" in v.reason  # keeps the message

    def test_provider_bug_surfaces(self):
        def buggy(_):
            raise AttributeError("typo in provider")  # code bug, not a runtime failure

        with pytest.raises(AttributeError):
            CofoldCrossCheck(provider=buggy).verify(_design(_site()))
