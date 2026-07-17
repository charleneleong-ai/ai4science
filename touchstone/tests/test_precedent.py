import math

import numpy as np
import pytest

from touchstone import BinderDesign, PrecedentHits, PrecedentVerifier, metalpdb_precedent_search
from touchstone.geometry.precedent import (
    MIN_MOTIF_HITS,
    MotifSelectivityVerifier,
    motif_enrichment,
)
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


class TestMotifSelectivity:
    """Metal discrimination from *observed occupancy* — the signal the geometry tiers can't give and
    the MLIP tier can't be trusted to give.

    Two things have to hold: it must recover metal/donor-set pairings biochemistry already settled,
    and it must not name a metal it has never actually observed on the donor set."""

    @pytest.mark.parametrize(
        "donors, expected, protein",
        [
            (("N", "N", "S", "S"), "Cu2+", "type-1 blue copper (His2/Cys/Met) — soft donors"),
            (("O",) * 6, "Mn2+", "all-oxygen — the hard ion wins, and Cu2+ loses"),
            (("S", "S", "S", "S"), "Fe3+", "Cys4 — iron-sulfur clusters dominate"),
        ],
    )
    def test_enrichment_recovers_known_metalloprotein_chemistry(self, donors, expected, protein):
        panel = ("Ni2+", "Cu2+", "Co2+", "Zn2+", "Fe3+", "Mn2+")
        enr = motif_enrichment(design(donors).site, panel)
        assert enr.preferred == expected, f"{protein}: got {enr.preferred}"

    def test_normalising_by_metal_abundance_is_what_makes_it_right(self):
        # the crux. Zn has ~10x more PDB sites than Cu, so RAW hits say Zn owns the type-1 copper
        # donor set (151 vs 17). Only after dividing by each metal's total presence does Cu win.
        panel = ("Cu2+", "Zn2+")
        enr = motif_enrichment(design(("N", "N", "S", "S")).site, panel)
        assert enr.hits["Zn2+"] > enr.hits["Cu2+"]           # raw counts point the wrong way ...
        assert enr.enrichment["Cu2+"] > enr.enrichment["Zn2+"]  # ... enrichment corrects them
        assert enr.preferred == "Cu2+"

    def test_trust_needs_the_target_to_actually_own_the_donor_set(self):
        # a Cu design on the soft-donor set trusts; the same site targeting Co does not
        soft = ("N", "N", "S", "S")
        v = MotifSelectivityVerifier(("Ni2+", "Cu2+", "Co2+"))
        assert v.verify(design(soft, metal="Cu2+")).trust
        assert not v.verify(design(soft, metal="Co2+")).trust

    def test_a_metal_never_seen_on_the_donor_set_cannot_own_it(self):
        # the smoothing trap: each metal's denominator is its own total, and those span 60x (Pt 53 vs
        # Zn 3105). A constant pseudo-count is therefore a *metal-dependent floor* — it handed rare
        # metals with ZERO hits a higher rate than abundant metals with real ones (25 of 85 donor
        # sets named a 0-hit metal). Shrinking toward the pooled background rate fixes it.
        panel = ("Ni2+", "Cu2+", "Co2+", "Zn2+", "Fe3+", "Mn2+", "Pt2+")
        for donors in [("O", "S", "S", "S", "S"), ("N", "N", "S", "S"), ("O",) * 6, ("S",) * 4]:
            enr = motif_enrichment(design(donors).site, panel)
            if enr.preferred is not None:
                assert enr.hits[enr.preferred] >= MIN_MOTIF_HITS, (
                    f"{enr.donors}: named {enr.preferred} the owner on {enr.hits[enr.preferred]} hits"
                )

    def test_defers_when_no_metal_has_enough_precedent(self):
        v = MotifSelectivityVerifier().verify(design(("S",) * 6))  # M-S6: not in the PDB
        assert v.ood and "too little evidence" in v.reason

    def test_a_thin_competitor_is_dropped_not_deferred(self):
        # Au3+ has 4 precedented sites. That's a hole in the reference table, not a fault in the
        # design — dropping Au from the panel beats deferring a design the user never asked about it.
        v = MotifSelectivityVerifier(("Ni2+", "Au3+")).verify(design(("N", "N", "S", "S"), metal="Cu2+"))
        assert not v.ood and v.trust  # still judged, on the Cu-vs-Ni comparison that survives
        assert "Au3+" in v.reason and "too few precedents" in v.reason

    def test_defers_when_the_target_itself_is_too_thin(self):
        v = MotifSelectivityVerifier(("Ni2+",)).verify(design(("N", "N", "S", "S"), metal="Au3+"))
        assert v.ood and "Au3+" in v.reason

    def test_a_panel_with_no_competitor_defers_rather_than_scoring_nan(self):
        # target-only panel: ratio has no competitor, and inf/(1+inf) is NaN — which also serialises
        # to invalid JSON and breaks --json for any non-Python consumer
        v = MotifSelectivityVerifier(("Ni2+",)).verify(design(("N", "N", "S", "S"), metal="Ni2+"))
        assert v.ood and "no competitor" in v.reason
        assert not math.isnan(v.score)

    def test_oxidation_states_of_one_element_are_not_competitors(self):
        # the table keys on the element, so Fe2+ and Fe3+ share counts — a panel holding both would
        # compare the target against itself (ratio 1.0, never trusting)
        v = MotifSelectivityVerifier(("Fe2+", "Zn2+"))
        panel = v.panel(design(("S",) * 4, metal="Fe3+"))
        assert panel == ["Fe3+", "Zn2+"]

    def test_zinc_is_a_real_competitor_for_the_type_1_copper_donor_set(self):
        # honest, and useful design guidance: N2S2 is only ~1.1x more Cu-characteristic than
        # Zn-characteristic (Cu 5.5%, Zn 4.9% — zinc fingers use Cys/His too), so a type-1 site is
        # NOT decisively Cu-selective against Zn. Cu still wins it against Ni/Co by ~3.8x.
        soft = ("N", "N", "S", "S")
        with_zn = MotifSelectivityVerifier(("Ni2+", "Co2+", "Zn2+")).verify(design(soft, metal="Cu2+"))
        without = MotifSelectivityVerifier(("Ni2+", "Co2+")).verify(design(soft, metal="Cu2+"))
        assert with_zn.metrics["preferred"] == "Cu2+" and not with_zn.trust  # preferred, not decisive
        assert without.trust
