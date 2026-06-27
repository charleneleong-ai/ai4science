import numpy as np
import pytest

from touchstone import BinderDesign, MogulVerifier
from touchstone.core import CoordinationSite
from touchstone.geometry.mogul import MogulFragment


def _design() -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.array([[2.1, 0, 0], [0, 2.1, 0]]), ("N", "O"))
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


def _analyse(z_scores, nhits=200):
    """A stand-in Mogul analyser returning fragments with the given z-scores."""

    def run(_site, _metal):
        return [MogulFragment(f"Ni-{i}", 2.1, z, nhits) for i, z in enumerate(z_scores)]

    return run


class TestMogulVerifier:
    def test_normal_geometry_trusts(self):
        v = MogulVerifier(_analyse([0.4, 0.9])).verify(_design())
        assert v.trust and not v.ood and v.score > 0.8

    def test_mildly_unusual_is_weak(self):
        # between trust_z and ood_z: judgeable, not trusted, not off-distribution
        v = MogulVerifier(_analyse([3.0])).verify(_design())
        assert not v.trust and not v.ood

    def test_off_distribution_bond_defers(self):
        v = MogulVerifier(_analyse([6.0])).verify(_design())  # past ood_z
        assert v.ood and not v.trust

    def test_thin_csd_support_defers(self):
        v = MogulVerifier(_analyse([0.5], nhits=3)).verify(_design())
        assert v.ood and "insufficient CSD support" in v.reason

    def test_no_fragments_defers(self):
        v = MogulVerifier(_analyse([])).verify(_design())
        assert v.ood and not v.trust

    def test_analyser_failure_defers(self):
        def boom(_site, _metal):
            raise RuntimeError("no CSD licence")

        v = MogulVerifier(boom).verify(_design())
        assert v.ood and "unavailable" in v.reason

    def test_worst_bond_drives_the_verdict(self):
        # one bad bond among good ones still fails trust (max |z|, not mean)
        v = MogulVerifier(_analyse([0.3, 0.5, 4.0])).verify(_design())
        assert not v.trust
