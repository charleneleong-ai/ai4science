import numpy as np
import pytest

from touchstone import BinderDesign, PrecedentHits, PrecedentVerifier, metalpdb_precedent_search
from touchstone.core import CoordinationSite

OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def design(elems=("N", "N", "O", "O", "N", "O"), metal: str = "Ni2+") -> BinderDesign:
    site = CoordinationSite(metal, np.zeros(3), OCT[: len(elems)] * 2.1, tuple(elems))
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


def searcher(nhits: int):
    def search(site) -> PrecedentHits:
        return PrecedentHits(nhits=nhits, motif="Ni-N3O3")

    return search


class TestPrecedentVerifier:
    def test_well_precedented_motif_trusts(self):
        v = PrecedentVerifier(searcher(40), min_hits=5).verify(design())
        assert v.trust and v.score > 0.5

    def test_thinly_precedented_is_weak(self):
        # 0 < hits < min_hits ⇒ judgeable but not trusted
        v = PrecedentVerifier(searcher(2), min_hits=5).verify(design())
        assert not v.trust and not v.ood

    def test_unprecedented_motif_defers(self):
        v = PrecedentVerifier(searcher(0)).verify(design())
        assert v.ood and "unprecedented" in v.reason

    def test_search_failure_defers(self):
        def boom(_site):
            raise RuntimeError("no licence")

        v = PrecedentVerifier(boom).verify(design())
        assert v.ood and "failed" in v.reason


class TestMetalPDBPrecedentSearch:
    """The open, licence-free default searcher — reads the bundled MetalPDB motif→count table."""

    def test_known_motif_has_precedents(self):
        hits = metalpdb_precedent_search(design(("N", "N", "N", "O", "O", "O")).site)  # Ni-N3O3
        assert hits.motif == "Ni-N3O3" and hits.nhits > 0

    def test_unseen_motif_has_zero_precedents(self):
        hits = metalpdb_precedent_search(design(("S", "S", "S", "S", "S", "S")).site)  # Ni-S6: not in the PDB
        assert hits.nhits == 0

    def test_default_verifier_runs_the_open_search(self):
        # PrecedentVerifier() with no injected searcher judges against the bundled table
        v = PrecedentVerifier().verify(design(("N", "N", "N", "O", "O", "O")))  # Ni-N3O3, well-precedented
        assert v.trust and not v.ood

    @pytest.mark.parametrize(
        "metal,donors,protein",
        [
            ("Ni2+", ("S", "S", "S", "S"), "NiFe hydrogenase (Cys4)"),
            ("Zn2+", ("S", "S", "S", "S"), "zinc finger (Cys4)"),
            ("Zn2+", ("N", "N", "N", "O"), "carbonic anhydrase (His3)"),
        ],
    )
    def test_canonical_metalloprotein_motifs_are_precedented(self, metal, donors, protein):
        # the guard on the table's domain. Built *without* excluding solvent, the counts drift to
        # aqua-padded pseudo-motifs (Ni-O6 — hexaaqua nickel — scored 22) and starve the real ones:
        # Ni-S4 fell to 3 hits, so touchstone called the NiFe-hydrogenase site unprecedented.
        v = PrecedentVerifier().verify(design(donors, metal=metal))
        assert v.trust, f"{protein} scored as unprecedented: {v.reason}"
