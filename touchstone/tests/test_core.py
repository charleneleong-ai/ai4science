import numpy as np
import pytest

from touchstone import octahedral_site
from touchstone.core import BinderDesign, Verdict, oxidation_state


class TestOxidationState:
    @pytest.mark.parametrize(
        "label, expected",
        [("Ni2+", 2), ("Fe3+", 3), ("Cu+", 1), ("Na+", 1), ("Cl-", -1), ("O2-", -2)],
    )
    def test_parses_charge(self, label, expected):
        # bare-sign labels like 'Cu+' (cuprous) must read as ±1, not raise
        assert oxidation_state(label) == expected


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
