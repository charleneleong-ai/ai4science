"""The Kaggle notebook vendors the package logic because its sandbox has no internet.

These tests pin the copy against the original. Without them the two drift, the notebook
quietly scores differently from everything measured locally, and the local numbers stop
describing the submission.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from casmi.adducts import ADDUCTS as PKG_ADDUCTS
from casmi.adducts import neutral_mass as pkg_neutral_mass
from casmi.formula import MONOISOTOPIC as PKG_MONO
from casmi.formula import formula_mass as pkg_formula_mass
from casmi.fragments import cached_fragments as pkg_fragments

NOTEBOOK = Path(__file__).parent.parent / "notebooks" / "kaggle_submission.py"

spec = importlib.util.spec_from_file_location("kaggle_submission", NOTEBOOK)
nb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nb)

MOLECULES = [
    "CCCO", "CC(C)O", "c1ccccc1", "OCC1OC(O)C(O)C(O)C1O",
    "CCOC(C)=O", "CN(Cc1ncccc1C(=O)O)CC1(c2ccc(Br)cc2)CC1",
]


class TestConstantsMatch:
    def test_element_masses(self):
        assert nb.MONOISOTOPIC == PKG_MONO

    def test_adduct_table(self):
        assert nb.ADDUCTS == PKG_ADDUCTS

    @pytest.mark.parametrize("name", ["H_MASS", "MAX_BREAKS", "MAX_BONDS", "H_SHIFTS"])
    def test_fragment_tunables(self, name):
        # a divergent break budget or H tolerance silently changes what the notebook submits
        from casmi import fragments as pkg_frag

        assert getattr(nb, name) == getattr(pkg_frag, name)

    def test_proton_mass(self):
        from casmi.adducts import PROTON

        assert nb.PROTON == PROTON


class TestBehaviourMatches:
    # The two implementations sum element masses in different orders, so they differ in the
    # last bits (~1e-11). Matching uses a 0.01 Da tolerance and masses are rounded to 1e-6
    # before set membership, so 1e-9 is seven orders below anything that can change a
    # decision — asserting bit-identity would pin float association, not behaviour.
    MASS_TOL = 1e-9

    @pytest.mark.parametrize("formula", ["H2O", "C6H12O6", "C21H23BrN2O2", "C2H6O", "", "Xx3"])
    def test_formula_mass(self, formula):
        got, want = nb.formula_mass(formula), pkg_formula_mass(formula)
        if want is None:
            assert got is None
        else:
            assert got == pytest.approx(want, abs=self.MASS_TOL)

    @pytest.mark.parametrize("adduct", sorted(PKG_ADDUCTS))
    def test_neutral_mass(self, adduct):
        assert nb.neutral_mass(300.0, adduct) == pytest.approx(
            pkg_neutral_mass(300.0, adduct), abs=self.MASS_TOL
        )

    def test_unknown_adduct_is_none_in_both(self):
        assert nb.neutral_mass(300.0, "[M+made-up]+") is pkg_neutral_mass(300.0, "[M+made-up]+") is None

    @pytest.mark.parametrize("smiles", MOLECULES)
    def test_fragment_masses_identical(self, smiles):
        # the ranker's entire signal — a divergence here means the notebook ranks differently.
        # Both round to 1e-6 before setting, so exact comparison is right at this layer.
        assert set(nb.fragment_masses(smiles)) == set(pkg_fragments(smiles))

    @pytest.mark.parametrize("smiles", MOLECULES)
    def test_fragment_masses_are_sorted(self, smiles):
        # the notebook's any_within bisects, so an unsorted tuple would silently miss matches
        frags = nb.fragment_masses(smiles)
        assert list(frags) == sorted(frags)

    def test_ranking_order_matches_including_ties(self):
        """The ranking order, not just the scores.

        The first version of the notebook sorted a (score, mass_error, smiles) tuple with a
        bare .sort(), so ties fell through to the SMILES string and broke alphabetically.
        With ~29 isomers sharing one formula mass, ties dominate — 202/400 submission rows
        diverged from the package. Scores agreeing is not enough; the order is the output.
        """
        from casmi.fragments import rank_by_fragments
        from casmi.harness import Query, Structure

        # three candidates that explain nothing, so every score ties and only the tie-break
        # separates them; SMILES chosen so alphabetical order differs from mass order
        cands = [
            Structure("K1", "CCCO", 60.057515),
            Structure("K2", "CC(C)O", 60.057515),
            Structure("K3", "COC", 60.057515),
        ]
        query = Query("Q", "", "[M+H]+", 60.057515 + nb.PROTON)
        query.spectra.append(([9999.0], [1.0]))

        pkg_order = [s.smiles for s in rank_by_fragments(query, cands)]

        positive = True
        scored = []
        for s in cands:
            frags = nb.fragment_masses(s.smiles)
            best = max((nb.explained_fraction(m, i, frags, positive) for m, i in query.spectra),
                       default=0.0)
            scored.append((-best, abs(s.mass - (query.mass or 0.0)), s.smiles))
        scored.sort(key=lambda row: (row[0], row[1]))
        nb_order = [smi for _, _, smi in scored]

        assert nb_order == pkg_order

    def test_explained_fraction_matches(self):
        from casmi.fragments import explained_fraction as pkg_explained

        frags = nb.fragment_masses("CCCO")
        mzs = [m + nb.PROTON for m in frags[:3]] + [999.0]
        intens = [0.4, 0.3, 0.2, 0.1]
        assert nb.explained_fraction(mzs, intens, frags, True) == pytest.approx(
            pkg_explained(mzs, intens, frozenset(frags), True)
        )


class TestPopularityPrior:
    """Shortlisting and blending against the PubChem popularity prior."""

    def test_top_distinct_keeps_most_popular_per_connectivity(self):
        pop = np.array([1.0, 5.0, 3.0, 4.0, 2.0], dtype=np.float32)
        keys = np.array([b"A", b"A", b"B", b"", b""], dtype="S14")
        # A's best is row 1; keyless rows have no prior and are dropped
        assert nb.top_distinct(pop, keys, 1) == [1]
        assert nb.top_distinct(pop, keys, 10) == [1, 2]

    def test_blend_lets_popularity_outrank_a_small_spectral_edge(self):
        explained = np.array([0.50, 0.49, 0.10])
        pop = np.array([0.0, 12.0, 24.0])
        assert list(np.argsort(-nb.blend(explained, pop))) == [2, 1, 0]

    def test_blend_without_prior_preserves_spectral_order(self):
        explained = np.array([0.2, 0.7, 0.4, 0.4])
        order = np.argsort(-nb.blend(explained, np.zeros(4)), kind="stable")
        assert list(order) == [1, 2, 3, 0]

    def test_blend_constant_spectral_scores_is_zero(self):
        assert not nb.blend(np.full(3, 0.3), np.zeros(3)).any()

    def test_train_candidate_takes_popularity_from_outside_the_shortlist(self, monkeypatch):
        monkeypatch.setattr(nb, "SHORTLIST", 1)
        smiles = [b"CCO", b"OCC", b"COC"]
        tier = object.__new__(nb.PubChemTier)
        tier.mass = np.array([46.0, 46.0, 46.0])
        tier.off = np.cumsum([0] + [len(x) for x in smiles])
        tier.smiles = np.frombuffer(b"".join(smiles), dtype=np.uint8)
        tier.lsid = np.array([9.0, 2.0, 1.0], dtype=np.float16)
        tier.lpmid = np.zeros(3, dtype=np.float16)
        tier.ik14 = np.array([b"POPULAR", b"TRAINKEY", b"OTHER"], dtype="S14")
        entries = [("TRAINKEY", ("OCC", 46.0)), ("NOTINPUBCHEM", ("CCC", 46.0))]

        cands = nb.candidates_for(46.0, entries, [46.0, 46.0], tier)

        assert {k: c.pop for k, c in cands.items()} == {
            "POPULAR": 9.0, "TRAINKEY": 2.0, "NOTINPUBCHEM": 0.0}
