import numpy as np
import pytest

from touchstone import octahedral_site, provider_from
from touchstone.core import BinderDesign, Verdict


class TestCoordinationSite:
    def test_bond_lengths_match_construction(self):
        site = octahedral_site("Ni2+", bond=2.1)
        assert np.allclose(site.bond_lengths(), 2.1)

    def test_coordination_number(self):
        assert octahedral_site("Ni2+").coordination_number == 6

    def test_octahedral_angles_are_90_or_180(self):
        angles = np.sort(octahedral_site("Ni2+").bond_angles())
        # 12 cis pairs at 90°, 3 trans pairs at 180°
        assert np.allclose(angles[:12], 90.0, atol=1e-6)
        assert np.allclose(angles[12:], 180.0, atol=1e-6)

    @pytest.mark.parametrize("metal", ["Ni2+", "Cu2+"])
    def test_target_metal_reads_from_site(self, metal):
        d = BinderDesign("SEQ", octahedral_site(metal), generator="x", generator_confidence=0.5)
        assert d.target_metal == metal


class TestVerdictLabel:
    @pytest.mark.parametrize(
        "trust, ood, expected",
        [(True, False, "trust"), (False, False, "weak"), (False, True, "defer"), (True, True, "defer")],
    )
    def test_label(self, trust, ood, expected):
        # ood takes precedence — a site off the manifold defers even if otherwise "trusted"
        assert Verdict(0.5, trust=trust, ood=ood, reason="x").label == expected

    def test_defer_factory(self):
        v = Verdict.defer("no data")
        assert not v.trust and v.ood and v.score == 0.0 and v.reason == "no data — defer"
        assert Verdict.defer("weak signal", score=0.3).score == 0.3  # ood verdict can still carry a score


class TestProviderFrom:
    """The shared factory behind the expression / thermostability / co-fold providers."""

    def _design(self, source: str | None = None) -> BinderDesign:
        return BinderDesign("SEQ", octahedral_site("Ni2+"), generator="x", generator_confidence=0.5, source=source)

    def test_looks_up_by_sequence_and_transforms_present_values(self):
        assert provider_from({"SEQ": 3.0}, transform=lambda v: v * 2)(self._design()) == 6.0

    def test_missing_returns_none_without_transforming(self):
        assert provider_from({}, transform=lambda v: v * 2)(self._design()) is None

    def test_alternate_key(self):
        assert provider_from({"p.pdb": "hit"}, key="source")(self._design(source="p.pdb")) == "hit"
