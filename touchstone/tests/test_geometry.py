import numpy as np
import pytest

from touchstone import (
    BinderDesign,
    GeometryVerifier,
    MockGenerator,
    MockReference,
    octahedral_site,
    rank,
    under_leachate,
)
from touchstone.core import CoordinationSite


def _design(site: CoordinationSite, conf: float = 0.7) -> BinderDesign:
    return BinderDesign("SEQ", site, generator="test", generator_confidence=conf)


@pytest.fixture
def verifier() -> GeometryVerifier:
    return GeometryVerifier(MockReference())


# A near-ideal Ni2+ site vs. two ways to be wrong: stretched bonds, and wrong coordination number.
GOOD = _design(octahedral_site("Ni2+", bond=2.09))
STRAINED = _design(octahedral_site("Ni2+", bond=2.6))
WRONG_CN = _design(
    CoordinationSite("Ni2+", np.zeros(3), octahedral_site("Ni2+").ligand_xyz[:4], ("N", "N", "O", "O"))
)


class TestGeometryVerifier:
    def test_good_site_is_trusted(self, verifier):
        v = verifier.verify(GOOD)
        assert v.trust and not v.ood

    @pytest.mark.parametrize("bad", [STRAINED, WRONG_CN], ids=["strained", "wrong_cn"])
    def test_bad_site_ranks_below_good(self, verifier, bad):
        assert verifier.verify(GOOD).score > verifier.verify(bad).score
        assert not verifier.verify(bad).trust

    def test_rank_orders_by_score(self, verifier):
        ordered = rank([STRAINED, GOOD, WRONG_CN], verifier)
        assert ordered[0][0] is GOOD
        assert [v.score for _, v in ordered] == sorted((v.score for _, v in ordered), reverse=True)


class TestExtremeConditionOOD:
    """The angle: a site trusted in benign conditions goes OOD under extreme leachate."""

    def test_ood_crossover_under_leachate(self, verifier):
        assert not verifier.verify(GOOD).ood  # in-distribution
        assert verifier.verify(under_leachate(GOOD, bond_stretch=0.6)).ood  # pushed off-manifold

    def test_leachate_collapses_trust(self, verifier):
        assert verifier.verify(GOOD).trust
        assert not verifier.verify(under_leachate(GOOD, bond_stretch=0.6)).trust


class TestGeneratorBlindness:
    """Swapping the generator must not change the verifier's API or behaviour."""

    def test_mock_generator_output_verifies(self, verifier):
        designs = MockGenerator(seed=1).design("Ni2+", n=5)
        verdicts = [verifier.verify(d) for d in designs]
        assert all(v.trust for v in verdicts)  # near-ideal sites ⇒ all trusted

    def test_verifier_ignores_generator_confidence(self, verifier):
        site = octahedral_site("Ni2+")
        low = verifier.verify(_design(site, conf=0.01))
        high = verifier.verify(_design(site, conf=0.99))
        assert (low.score, low.trust, low.ood) == (high.score, high.trust, high.ood)


class TestReferencePluggability:
    def test_swapping_reference_changes_verdict_not_api(self):
        """A reference whose ideal bond is far from the site flips trust — same call."""
        loose = MockReference()
        loose._TABLE = {"Ni2+": loose._TABLE["Ni2+"].__class__("Ni2+", 6, 3.5, 0.05)}
        assert GeometryVerifier(MockReference()).verify(GOOD).trust
        assert not GeometryVerifier(loose).verify(GOOD).trust

    def test_unknown_metal_raises(self, verifier):
        with pytest.raises(KeyError):
            verifier.verify(_design(octahedral_site("Au3+")))


class TestEmptySite:
    def test_no_donors_defers_with_finite_score(self, verifier):
        """A site with no coordinating atoms (CN=0) is the worst case, not NaN."""
        empty = _design(CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ()))
        v = verifier.verify(empty)
        assert v.score == 0.0 and not v.trust and v.ood
        assert rank([empty, GOOD], verifier)[0][0] is GOOD  # NaN would break this sort
