from touchstone import (
    BinderDesign,
    GeometryVerifier,
    MockReference,
    octahedral_site,
    stress_profile,
    under_low_pH,
)

V = GeometryVerifier(MockReference())


def _design(site) -> BinderDesign:
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


class TestStressProfile:
    def test_returns_the_condition_map(self):
        prof = stress_profile(_design(octahedral_site("Ni2+", bond=2.09)), V)
        assert set(prof) == {"neutral", "leachate", "low_pH"}

    def test_ideal_site_trusts_at_baseline_but_not_under_leachate(self):
        prof = stress_profile(_design(octahedral_site("Ni2+", bond=2.09)), V)
        assert prof["neutral"].trust and not prof["leachate"].trust  # stretched bonds ⇒ off-manifold

    def test_severe_protonation_collapses_trust(self):
        # losing one donor from CN6 is survivable (the robustness point); losing three
        # drops below the observed CN range → no longer trusted
        prof = stress_profile(_design(octahedral_site("Ni2+", bond=2.09)), V, n_protonate=3)
        assert not prof["low_pH"].trust


def test_protonation_drops_the_n_weakest_donors():
    site = octahedral_site("Ni2+", bond=2.09)  # CN 6
    assert under_low_pH(_design(site), n_drop=2).site.coordination_number == 4
