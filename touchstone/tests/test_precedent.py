import numpy as np
import pytest

from touchstone import BinderDesign, PrecedentHits, PrecedentVerifier
from touchstone.core import CoordinationSite

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def _design(elems=("N", "N", "O", "O", "N", "O")) -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), _OCT[: len(elems)] * 2.1, tuple(elems))
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


def _searcher(nhits: int):
    def search(site) -> PrecedentHits:
        return PrecedentHits(nhits=nhits, motif="Ni-N3O3")

    return search


class TestPrecedentVerifier:
    def test_well_precedented_motif_trusts(self):
        v = PrecedentVerifier(_searcher(40), min_hits=5).verify(_design())
        assert v.trust and v.score > 0.5

    def test_thinly_precedented_is_weak(self):
        # 0 < hits < min_hits ⇒ judgeable but not trusted
        v = PrecedentVerifier(_searcher(2), min_hits=5).verify(_design())
        assert not v.trust and not v.ood

    def test_unprecedented_motif_defers(self):
        v = PrecedentVerifier(_searcher(0)).verify(_design())
        assert v.ood and "unprecedented" in v.reason

    def test_search_failure_defers(self):
        def boom(_site):
            raise RuntimeError("no licence")

        v = PrecedentVerifier(boom).verify(_design())
        assert v.ood and "failed" in v.reason
