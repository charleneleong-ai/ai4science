import pytest

from casmi.adducts import PROTON
from casmi.fragments import (
    H_MASS,
    any_within,
    cached_fragments,
    explained_fraction,
    fragment_mass_set,
    rank_by_fragments,
    score_candidate,
)
from casmi.harness import Query, Structure

# C3H8O isomers that differ only by branching. Breaking analogous bonds gives pieces of the
# SAME formula, so their fragment-mass sets are identical — the method's blind spot.
PROPANOL = "CCCO"
ISOPROPANOL = "CC(C)O"

# C4H8O2 isomers that differ in where the ester oxygen sits, so cleavage gives genuinely
# different fragment formulas (acetate 59.013 + ethyl vs propanoate 73.029 + methyl).
ETHYL_ACETATE = "CCOC(C)=O"
METHYL_PROPANOATE = "CCC(=O)OC"

GLUCOSE = "OCC1OC(O)C(O)C(O)C1O"


class TestFragmentEnumeration:
    def test_intact_molecule_is_always_a_fragment(self):
        frags = fragment_mass_set(PROPANOL)
        assert any(abs(m - 60.057515) < 1e-3 for m in frags)  # C3H8O

    def test_breaking_bonds_yields_smaller_pieces(self):
        frags = fragment_mass_set(PROPANOL)
        assert min(frags) < 60.0, "no sub-fragment was produced"

    def test_isomers_give_different_fragment_sets(self):
        # the premise of this ranker: isomers that cleave into different formulas are separable
        assert fragment_mass_set(ETHYL_ACETATE) != fragment_mass_set(METHYL_PROPANOATE)

    @pytest.mark.parametrize("breaks", [1, 2])
    def test_branching_isomers_are_separable(self, breaks):
        # An earlier version reported these as identical, which was a hydrogen-accounting bug
        # in RDKit's fragmenter (addDummies=False + sanitizeFrags=False leaves stale implicit
        # H counts), not a property of fragmentation. Pinned so that regression is visible:
        # propanol keeps C2H5/CH3O at 29.039/31.018, which isopropanol cannot produce.
        assert fragment_mass_set(PROPANOL, max_breaks=breaks) != fragment_mass_set(ISOPROPANOL, max_breaks=breaks)

    def test_fragment_masses_are_deduplicated(self):
        # summing components in different orders differs in the last bits; a float set would
        # keep those as separate entries and bloat every lookup
        frags = fragment_mass_set("c1ccccc1", max_breaks=1)
        assert len(frags) == 1

    def test_unparseable_smiles_returns_empty_not_an_exception(self):
        assert fragment_mass_set("not-a-smiles") == frozenset()

    def test_more_breaks_never_loses_fragments(self):
        one = fragment_mass_set(GLUCOSE, max_breaks=1)
        two = fragment_mass_set(GLUCOSE, max_breaks=2)
        assert one <= two

    def test_ring_cleavage_finds_fragments_that_acyclic_only_cannot(self):
        # the fix: 68% of bonds in these natural products are ring bonds, and a ring needs
        # two cuts to open, so excluding them left most of the spectrum unexplainable
        acyclic_only = fragment_mass_set(GLUCOSE, include_rings=False)
        with_rings = fragment_mass_set(GLUCOSE, include_rings=True)
        assert acyclic_only < with_rings          # strictly more reachable, nothing lost
        assert len(with_rings) > 2 * len(acyclic_only)  # measured ~2.9x on glucose

    def test_one_ring_cut_alone_cannot_disconnect_a_ring(self):
        # benzene: a single cut leaves it in one piece, so only the intact mass appears
        frags = fragment_mass_set("c1ccccc1", max_breaks=1)
        assert len(frags) == 1
        assert next(iter(frags)) == pytest.approx(78.046950, abs=1e-4)

    def test_two_ring_cuts_open_an_aromatic_ring(self):
        # aromatic bonds must be breakable — the first version filtered to BondType.SINGLE
        # and so could never cut a benzene ring at all
        assert len(fragment_mass_set("c1ccccc1", max_breaks=2)) > 1

    def test_cache_returns_the_same_object(self):
        assert cached_fragments(PROPANOL) is cached_fragments(PROPANOL)


class TestMatching:
    def test_exact_match_is_found(self):
        assert any_within([100.0, 200.0], 100.0, 0.01)

    def test_outside_tolerance_is_rejected(self):
        assert not any_within([100.0, 200.0], 150.0, 0.01)

    def test_hydrogen_rearrangement_is_allowed(self):
        # fragments commonly appear +-1 H from the bare mass
        assert any_within([100.0], 100.0 + H_MASS, 0.01)
        assert any_within([100.0], 100.0 - H_MASS, 0.01)

    def test_two_hydrogens_off_is_not_matched(self):
        assert not any_within([100.0], 100.0 + 2 * H_MASS, 0.01)


class TestExplainedFraction:
    def test_all_intensity_explained_scores_one(self):
        frags = frozenset({100.0})
        assert explained_fraction([100.0 + PROTON], [1.0], frags, positive=True) == pytest.approx(1.0)

    def test_nothing_explained_scores_zero(self):
        assert explained_fraction([500.0], [1.0], frozenset({100.0}), positive=True) == 0.0

    def test_it_is_intensity_weighted_not_a_peak_count(self):
        # one explained base peak should outweigh an unexplained minor peak
        frags = frozenset({100.0})
        got = explained_fraction([100.0 + PROTON, 400.0], [0.9, 0.1], frags, positive=True)
        assert got == pytest.approx(0.9)

    def test_negative_mode_uses_the_other_charge(self):
        frags = frozenset({100.0})
        assert explained_fraction([100.0 - PROTON], [1.0], frags, positive=False) == pytest.approx(1.0)

    def test_empty_spectrum_scores_zero(self):
        assert explained_fraction([], [], frozenset({100.0}), positive=True) == 0.0


class TestRanking:
    def query_for(self, smiles_mass: float, peaks: list[float]) -> Query:
        q = Query("TRUTH", "", "[M+H]+", smiles_mass + PROTON)
        q.spectra.append((peaks, [1.0] * len(peaks)))
        return q

    def test_structure_explaining_the_peaks_ranks_first(self):
        # The peak must be one ethyl acetate cannot reach *even allowing +-1 H* — set
        # subtraction alone is not enough, since the H tolerance lets a near-miss match.
        ea = sorted(fragment_mass_set(ETHYL_ACETATE))
        only_mp = [m for m in sorted(fragment_mass_set(METHYL_PROPANOATE) - set(ea))
                   if not any_within(ea, m, 0.01)]
        assert only_mp, "fixture invalid: every MP fragment is reachable by EA"
        q = self.query_for(88.052429, [only_mp[0] + PROTON])
        cands = [Structure("EA", ETHYL_ACETATE, 88.052429),
                 Structure("MP", METHYL_PROPANOATE, 88.052429)]
        assert rank_by_fragments(q, cands)[0].inchikey14 == "MP"

    def test_ties_fall_back_to_mass_error_not_arbitrary_order(self):
        q = self.query_for(60.057515, [999.0])  # nothing explicable
        near = Structure("NEAR", PROPANOL, 60.057515)
        far = Structure("FAR", PROPANOL, 60.0590)
        assert rank_by_fragments(q, [far, near])[0].inchikey14 == "NEAR"

    def test_unparseable_candidate_scores_zero_and_does_not_crash(self):
        q = self.query_for(60.057515, [43.0 + PROTON])
        assert score_candidate(q, Structure("BAD", "not-a-smiles", 60.0)) == 0.0

    def test_score_takes_the_best_spectrum_not_the_mean(self):
        # a low-energy spectrum that barely fragments must not penalise a good structure
        q = Query("T", "", "[M+H]+", 60.057515 + PROTON)
        q.spectra.append(([999.0], [1.0]))                     # explains nothing
        q.spectra.append(([60.057515 + PROTON], [1.0]))        # explains everything
        assert score_candidate(q, Structure("P", PROPANOL, 60.057515)) == pytest.approx(1.0)
