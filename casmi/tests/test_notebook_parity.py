"""The Kaggle notebook vendors the package logic because its sandbox has no internet.

These tests pin the copy against the original. Without them the two drift, the notebook
quietly scores differently from everything measured locally, and the local numbers stop
describing the submission.
"""

import importlib.util
from collections.abc import Callable, Iterable
from pathlib import Path
from types import ModuleType
from typing import Any

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


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def trainset():
    return load_module("ranker_trainset", NOTEBOOK.parent / "ranker" / "ranker_trainset.py")


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

    def test_candidates_drop_excluded_keys(self):
        entries = [("PARENT", ("CCO", 46.0)), ("OTHER", ("OCC", 46.0))]
        cands = nb.candidates_for(46.0, entries, [46.0, 46.0], None, {"PARENT"})
        assert list(cands) == ["OTHER"]

    def test_copied_parents_match_bit_identical_spectra_only(self, tmp_path):
        import pyarrow as pa
        import pyarrow.parquet as pq

        rows = [  # (inchikey14, adduct, precursor, mzs, intensities)
            ("COPIED", "[M+H]+", 100.0, [50.0, 60.0], [1.0, 0.5]),
            ("NEARMISS", "[M+H]+", 100.0, [50.0, 60.0], [1.0, 0.4999]),
            ("OTHERADDUCT", "[M+Na]+", 100.0, [50.0, 60.0], [1.0, 0.5]),
        ]
        cols = ["inchikey14", "adduct", "precursor_mz", "ms2_mzs", "ms2_normalized_intensities"]
        pq.write_table(pa.table({c: [r[i] for r in rows] for i, c in enumerate(cols)}), tmp_path / "t.parquet")
        mol = {"records": [dict(mz=[50.0, 60.0], it=[1.0, 0.5], prec=100.0, adduct="[M+H]+")]}

        assert nb.copied_parents(tmp_path / "t.parquet", {"m1": mol}) == {"m1": {"COPIED"}}


class TestRankerTrainset:
    """Features and labels the ranker's training-set kernel builds per query."""

    def test_rows_label_the_query_and_rank_within_it(self, trainset, monkeypatch):
        monkeypatch.setattr(nb, "SHORTLIST", 3)
        smiles = [b"CCO", b"OCC", b"COC", b"CCCC"]
        tier = object.__new__(nb.PubChemTier)
        tier.mass = np.array([46.0, 46.0, 46.0, 58.0])
        tier.off = np.cumsum([0] + [len(x) for x in smiles])
        tier.smiles = np.frombuffer(b"".join(smiles), dtype=np.uint8)
        tier.lsid = np.array([1.0, 3.0, 2.0, 9.0], dtype=np.float16)
        tier.lpmid = np.zeros(4, dtype=np.float16)
        tier.ik14 = np.array([b"OTHER", b"QUERYKEY", b"THIRD", b"FAR"], dtype="S14")
        prior = object.__new__(trainset.Prior)
        prior.lpatent = np.zeros(4, dtype=np.float16)
        prior.np_bits = np.array([0, 5, 2, 0], dtype=np.uint8)

        class FakeFPNet:
            def scores(self, records, smis):
                return np.array([{"CCO": 0.1, "OCC": 0.9, "COC": 0.5}[s] for s in smis])

            def forget(self):
                pass

        mol = {"adduct": "[M+H]+", "records": [dict(mz=[30.0], it=[1.0], prec=47.0, adduct="[M+H]+",
                                                     mode=1, ce=None, ce_n=1)]}
        rows = trainset.candidate_rows(nb, tier, prior, "QUERYKEY", mol, 46.0, FakeFPNet())

        (query,) = [r for r in rows if r["label"]]
        assert (query["fz_rank"], query["pop_rank"]) == (0, 0) and query["fz_z"] > 0
        assert (query["np_coconut"], query["np_lotus"], query["np_npatlas"]) == (1, 0, 1)
        assert all(r["n_isomers"] == 3 and r["n_candidates"] == 3 for r in rows)  # 58 Da row is outside 5 ppm


class SerialPool:
    def map(self, fn: Callable[[Any], Any], items: Iterable[Any], chunksize: int | None = None) -> list[Any]:
        return [fn(x) for x in items]


class TestRecallProbe:
    """The second-pass probe reaches a structure the popularity shortlist drops."""

    def test_spectral_pass_reaches_an_unpopular_answer(self, monkeypatch):
        from rdkit import Chem
        from rdkit.Chem import rdFingerprintGenerator

        probe = load_module("recall_probe", NOTEBOOK.parent / "recall" / "recall_probe.py")
        monkeypatch.setattr(nb, "SHORTLIST", 1)
        smiles = [b"COC", b"CCO", b"CCCC"]
        tier = object.__new__(nb.PubChemTier)
        tier.mass = np.array([46.0, 46.0, 58.0])
        tier.off = np.cumsum([0] + [len(x) for x in smiles])
        tier.smiles = np.frombuffer(b"".join(smiles), dtype=np.uint8)
        tier.lsid = np.array([9.0, 1.0, 9.0], dtype=np.float16)
        tier.lpmid = np.zeros(3, dtype=np.float16)
        tier.ik14 = np.array([b"OTHER", b"QUERYKEY", b"FAR"], dtype="S14")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=probe.ECFP_SIZE)
        morgan = lambda s: gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))  # noqa: E731
        bits = np.flatnonzero(morgan("CCO") & ~morgan("COC"))
        chem = type("Chem", (), {"raw_fingerprint": staticmethod(morgan)})
        monkeypatch.setattr(probe, "W", dict(tier=tier, chem=chem, bits=bits, Chem=Chem, gen=gen,
                                            ecfp_cols=np.arange(len(bits)), ecfp_raw=bits))
        stats = {"truth_in_window": 0, "window_rows": []}

        rows = probe.probe_rows(nb, SerialPool(), "QUERYKEY", 46.0, np.ones(len(bits), np.float32), stats)

        by_key = {r["key"]: r for r in rows}
        assert set(by_key) == {"OTHER", "QUERYKEY"}  # 58 Da row is outside 5 ppm
        assert by_key["OTHER"]["in_pop"] and not by_key["QUERYKEY"]["in_pop"]
        assert by_key["QUERYKEY"]["label"] and by_key["QUERYKEY"]["pass1_rank"] == 0
        assert by_key["QUERYKEY"]["fz"] > by_key["OTHER"]["fz"]
        assert stats == {"truth_in_window": 1, "window_rows": [2]}

    def test_window_without_keyed_rows_is_a_counted_miss(self, monkeypatch):
        probe = load_module("recall_probe", NOTEBOOK.parent / "recall" / "recall_probe.py")
        tier = object.__new__(nb.PubChemTier)
        tier.mass, tier.ik14 = np.array([58.0]), np.array([b""], dtype="S14")
        monkeypatch.setattr(probe, "W", {"tier": tier})
        stats = {"truth_in_window": 0, "window_rows": []}
        for mass in (46.0, 58.0):  # empty window, then a window whose only row has no key
            assert probe.probe_rows(nb, SerialPool(), "QUERYKEY", mass, np.ones(1), stats) == []
        assert stats == {"truth_in_window": 0, "window_rows": [0, 1]}


class TestSecondPass:
    """The notebook's spectrum-driven pass adds what the shortlist missed, once per connectivity."""

    def test_adds_best_new_connectivities_only(self, monkeypatch):
        from rdkit.Chem import rdFingerprintGenerator

        smiles = [b"COC", b"CCO", b"OCC", b"CC=O", b"C=CO"]
        tier = object.__new__(nb.PubChemTier)
        tier.mass = np.array([46.0, 46.0, 46.0, 46.0, 46.0])
        tier.off = np.cumsum([0] + [len(x) for x in smiles])
        tier.smiles = np.frombuffer(b"".join(smiles), dtype=np.uint8)
        tier.lsid, tier.lpmid = np.zeros(5, dtype=np.float16), np.zeros(5, dtype=np.float16)
        tier.ik14 = np.array([b"TAKEN", b"ANSWER", b"ANSWER", b"", b"OTHER"], dtype="S14")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=nb.ECFP_SIZE)
        morgan = lambda s: gen.GetFingerprintAsNumPy(nb.Chem.MolFromSmiles(s))  # noqa: E731
        bits = np.flatnonzero(morgan("CCO"))
        monkeypatch.setattr(nb, "PASS_STATE", dict(tier=tier, bits=bits, ecfp_raw=bits, gen=gen))
        monkeypatch.setattr(nb, "pass_chem", lambda: type("Chem", (), {"raw_fingerprint": staticmethod(morgan)}))
        monkeypatch.setattr(nb, "SECOND_PASS_KEEP", 1)
        second = object.__new__(nb.SecondPass)
        second.tier, second.cols, second.nbits, second.pool = tier, np.arange(len(bits)), len(bits), SerialPool()

        added = second.candidates(45.9, 46.1, np.ones(len(bits), np.float32), {"TAKEN"})

        assert list(added) == ["ANSWER"] and added["ANSWER"].smiles in ("CCO", "OCC")
        assert second.candidates(50.0, 51.0, np.ones(len(bits), np.float32), set()) == {}
        monkeypatch.setattr(nb, "SECOND_PASS_N1", 1)  # cut before filtering: the best row is taken
        assert second.candidates(45.9, 46.1, np.ones(len(bits), np.float32), {"ANSWER"}) == {}
