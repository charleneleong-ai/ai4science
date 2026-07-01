import numpy as np
import pytest

from touchstone import BinderDesign, BondValenceVerifier
from touchstone.core import CoordinationSite

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def _design(metal: str = "Ni2+", bond: float = 2.06, n: int = 6, elem: str = "O") -> BinderDesign:
    site = CoordinationSite(metal, np.zeros(3), _OCT[:n] * bond, (elem,) * n)
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


@pytest.fixture
def verifier() -> BondValenceVerifier:
    return BondValenceVerifier()


class TestBondValenceVerifier:
    def test_charge_recovering_site_is_trusted(self, verifier):
        # 6 Ni–O at ~2.06 Å sum to a bond-valence of ≈2 = the formal charge of Ni²⁺
        v = verifier.verify(_design(bond=2.06))
        assert v.trust and not v.ood and v.score > 0.8

    @pytest.mark.parametrize("bond", [2.6, 1.8], ids=["under-bonded", "over-bonded"])
    def test_mischarged_site_defers(self, verifier, bond):
        # stretched bonds undershoot the valence sum; crushed bonds overshoot — both implausible
        v = verifier.verify(_design(bond=bond))
        assert not v.trust and v.ood

    def test_unknown_metal_raises(self, verifier):
        with pytest.raises(KeyError):
            verifier.verify(_design(metal="Au3+"))

    def test_empty_site_defers_with_finite_score(self, verifier):
        empty = BinderDesign(
            "SEQ", CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ()),
            generator="test", generator_confidence=0.5,
        )
        v = verifier.verify(empty)
        assert v.score == 0.0 and v.ood and not v.trust

    def test_unparameterized_donor_defers(self, verifier):
        # P is not in the R0 table — must defer, not silently undercount the BVS
        v = verifier.verify(_design(elem="P"))
        assert v.ood and not v.trust and "not parameterized" in v.reason
