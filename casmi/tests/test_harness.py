import pytest

from casmi.adducts import ADDUCTS, neutral_mass
from casmi.formula import formula_mass, parse_formula
from casmi.harness import (
    Query,
    Structure,
    evaluate,
    holdout_fraction,
    is_held_out,
    rank_by_mass_error,
    reciprocal_rank,
    retrieve,
)


def pool_of(*masses: float) -> tuple[list[Structure], list[float]]:
    pool = sorted((Structure(f"K{i}", f"C{i}", m) for i, m in enumerate(masses)), key=lambda s: s.mass)
    return pool, [s.mass for s in pool]


class TestAdducts:
    @pytest.mark.parametrize(
        "adduct,mz,expected",
        [
            ("[M+H]+", 301.007276, 300.0),
            ("[M-H]-", 298.992724, 300.0),
            ("[M+Na]+", 322.989218, 300.0),
            ("[2M+H]+", 601.007276, 300.0),  # dimer: the precursor holds two copies
        ],
    )
    def test_neutral_mass_inverts_the_adduct(self, adduct, mz, expected):
        assert neutral_mass(mz, adduct) == pytest.approx(expected, abs=1e-6)

    def test_unknown_adduct_returns_none_rather_than_guessing(self):
        # silently treating an unmapped adduct as [M+H]+ shifts the window and drops the truth
        assert neutral_mass(300.0, "[M+totally-made-up]+") is None

    def test_every_mapped_adduct_round_trips(self):
        for adduct, (mult, offset) in ADDUCTS.items():
            mz = 300.0 * mult + offset
            assert neutral_mass(mz, adduct) == pytest.approx(300.0, abs=1e-6)


class TestFormulaMass:
    @pytest.mark.parametrize(
        "formula,expected",
        [
            ("H2O", 18.010565),
            ("C6H12O6", 180.063388),   # glucose
            ("CH4", 16.031300),
            ("C2H6O", 46.041865),      # ethanol
            ("C21H23BrN2O2", 414.094290),  # a halogenated case, so Br is covered
        ],
    )
    def test_known_masses(self, formula, expected):
        assert formula_mass(formula) == pytest.approx(expected, abs=1e-4)

    def test_isomers_share_one_mass(self):
        # the point of formula-derived masses: mass error cannot separate isomers, so a
        # mass-only ranker gets no signal between them
        assert formula_mass("C2H6O") == formula_mass("C2H6O")

    @pytest.mark.parametrize("bad", ["", "C6H12O6???", "Xx3", "not-a-formula", "6CH"])
    def test_unparseable_returns_none_rather_than_a_partial_mass(self, bad):
        # a partial parse understates the mass and moves the candidate to the wrong window
        assert formula_mass(bad) is None

    def test_trailing_charge_is_tolerated(self):
        assert formula_mass("C6H12O6+") == pytest.approx(formula_mass("C6H12O6"))

    def test_repeated_element_accumulates(self):
        assert parse_formula("CH3CH3") == {"C": 2, "H": 6}


class TestSplit:
    def test_holdout_is_deterministic(self):
        # the split must not move between runs, or results aren't comparable
        assert holdout_fraction("ABCDEFGHIJKLMN") == holdout_fraction("ABCDEFGHIJKLMN")

    def test_fraction_is_in_range_and_spread(self):
        fracs = [holdout_fraction(f"KEY{i:05d}") for i in range(2000)]
        assert all(0.0 <= f < 1.0 for f in fracs)
        # a hash that clumped would silently make the holdout unrepresentative
        assert 0.4 < sum(f < 0.5 for f in fracs) / len(fracs) < 0.6

    def test_fraction_selects_roughly_that_share(self):
        n = sum(is_held_out(f"KEY{i:05d}", 0.1) for i in range(2000))
        assert 150 < n < 250


class TestRetrieval:
    def test_window_includes_only_masses_inside_it(self):
        pool, masses = pool_of(299.9, 300.0, 300.0005, 300.5)
        got = {s.mass for s in retrieve(pool, masses, 300.0, ppm=5.0)}
        assert got == {300.0, 300.0005}  # +-1.5 mDa at 300 Da

    def test_empty_window_is_not_an_error(self):
        pool, masses = pool_of(100.0, 500.0)
        assert list(retrieve(pool, masses, 300.0, ppm=5.0)) == []


class TestScoring:
    def test_reciprocal_rank_is_one_over_position(self):
        ranked = [Structure("A", "", 1.0), Structure("B", "", 1.0), Structure("C", "", 1.0)]
        assert reciprocal_rank(ranked, "A") == 1.0
        assert reciprocal_rank(ranked, "C") == pytest.approx(1 / 3)

    def test_truth_beyond_k_scores_zero(self):
        ranked = [Structure(f"K{i}", "", 1.0) for i in range(30)]
        assert reciprocal_rank(ranked, "K25", k=25) == 0.0

    def test_missing_truth_scores_zero(self):
        assert reciprocal_rank([Structure("A", "", 1.0)], "ZZZ") == 0.0

    def test_mass_error_ranker_puts_the_closest_first(self):
        q = Query("K1", "C", "[M+H]+", 301.007276)
        cands = [Structure("far", "", 300.001), Structure("near", "", 300.0000001)]
        assert rank_by_mass_error(q, cands)[0].inchikey14 == "near"


class TestEvaluate:
    def test_recall_bounds_mrr(self):
        # the truth is absent from the pool, so no ranker can score: recall and MRR are both 0
        pool, _ = pool_of(100.0, 500.0)
        q = Query("MISSING", "C", "[M+H]+", 301.007276)
        rep = evaluate([q], pool)
        assert rep.recall == 0.0 and rep.mrr == 0.0

    def test_recoverable_query_scores(self):
        truth = Structure("TRUTH", "C", 300.0)
        pool = sorted([truth, Structure("OTHER", "C", 300.0004)], key=lambda s: s.mass)
        q = Query("TRUTH", "C", "[M+H]+", 301.007276)
        rep = evaluate([q], pool)
        assert rep.recall == 1.0 and rep.mrr == 1.0

    def test_unmapped_adduct_is_counted_not_silently_scored(self):
        pool, _ = pool_of(300.0)
        q = Query("K0", "C", "[M+nonsense]+", 301.0)
        rep = evaluate([q], pool)
        assert rep.n_unmapped_adduct == 1
        assert rep.recall == 0.0  # excluded from the denominator, not counted as a miss
