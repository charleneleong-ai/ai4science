import json

import numpy as np
import pytest

from touchstone import (
    BinderDesign,
    CSDReference,
    GeometryVerifier,
    MetalHawkPrediction,
    MetalHawkVerifier,
    MetalPDBReference,
    MockGenerator,
    MockReference,
    PDBReference,
    octahedral_site,
    rank,
    selectivity_profile,
    under_leachate,
)
from touchstone.core import CoordinationSite
from touchstone.geometry import reference as _ref_mod
from touchstone.geometry.metalhawk import load_predictions, score_provider
from touchstone.geometry.reference import CSD_DATA, METALPDB_DATA, best_reference


def _boom(_design):  # a scorer that fails — module-level so it's not a lambda
    raise RuntimeError("metalhawk exploded")


def _design(site: CoordinationSite, conf: float = 0.7) -> BinderDesign:
    return BinderDesign("SEQ", site, generator="test", generator_confidence=conf)


@pytest.fixture
def verifier() -> GeometryVerifier:
    return GeometryVerifier(MockReference())


# A near-ideal Ni2+ site vs. two ways to be wrong: stretched bonds, and wrong coordination number.
GOOD = _design(octahedral_site("Ni2+", bond=2.09))
STRAINED = _design(octahedral_site("Ni2+", bond=2.6))
WRONG_CN = _design(
    CoordinationSite("Ni2+", np.zeros(3), octahedral_site("Ni2+").ligand_xyz[:4], ("N", "N", "O", "O"))
)


class TestGeometryVerifier:
    def test_good_site_is_trusted(self, verifier):
        v = verifier.verify(GOOD)
        assert v.trust and not v.ood

    @pytest.mark.parametrize("bad", [STRAINED, WRONG_CN], ids=["strained", "wrong_cn"])
    def test_bad_site_ranks_below_good(self, verifier, bad):
        assert verifier.verify(GOOD).score > verifier.verify(bad).score
        assert not verifier.verify(bad).trust

    def test_rank_orders_by_score(self, verifier):
        ordered = rank([STRAINED, GOOD, WRONG_CN], verifier)
        assert ordered[0][0] is GOOD
        assert [v.score for _, v in ordered] == sorted((v.score for _, v in ordered), reverse=True)


class TestExtremeConditionOOD:
    """The angle: a site trusted in benign conditions goes OOD under extreme leachate."""

    def test_ood_crossover_under_leachate(self, verifier):
        assert not verifier.verify(GOOD).ood  # in-distribution
        assert verifier.verify(under_leachate(GOOD, bond_stretch=0.6)).ood  # pushed off-manifold

    def test_leachate_collapses_trust(self, verifier):
        assert verifier.verify(GOOD).trust
        assert not verifier.verify(under_leachate(GOOD, bond_stretch=0.6)).trust


class TestGeneratorBlindness:
    """Swapping the generator must not change the verifier's API or behaviour."""

    def test_mock_generator_output_verifies(self, verifier):
        designs = MockGenerator(seed=1).design("Ni2+", n=5)
        verdicts = [verifier.verify(d) for d in designs]
        assert all(v.trust for v in verdicts)  # near-ideal sites ⇒ all trusted

    def test_verifier_ignores_generator_confidence(self, verifier):
        site = octahedral_site("Ni2+")
        low = verifier.verify(_design(site, conf=0.01))
        high = verifier.verify(_design(site, conf=0.99))
        assert (low.score, low.trust, low.ood) == (high.score, high.trust, high.ood)


class TestReferencePluggability:
    def test_swapping_reference_changes_verdict_not_api(self):
        """A reference whose ideal bond is far from the site flips trust — same call."""
        loose = MockReference()
        loose._TABLE = {"Ni2+": loose._TABLE["Ni2+"].__class__("Ni2+", 6, 3.5, 0.05)}
        assert GeometryVerifier(MockReference()).verify(GOOD).trust
        assert not GeometryVerifier(loose).verify(GOOD).trust

    def test_unknown_metal_raises(self, verifier):
        with pytest.raises(KeyError):
            verifier.verify(_design(octahedral_site("Au3+")))


class TestPDBReference:
    """The real reference, loaded from the cached PDB pull."""

    @pytest.mark.parametrize("metal", ["Ni2+", "Cu2+", "Co2+"])
    def test_geometry_is_physically_sane(self, metal):
        g = PDBReference().geometry(metal)
        assert 1.8 <= g.bond_length_mean <= 2.6  # metal–donor bonds live here
        assert g.bond_length_std > 0 and g.coordination_number >= 3

    def test_unknown_metal_raises(self):
        with pytest.raises(KeyError):
            PDBReference().geometry("Au3+")

    def test_clean_site_at_empirical_geometry_trusts(self):
        """A site matching the real PDB geometry (modal CN, mean bond) trusts."""
        g = PDBReference().geometry("Ni2+")
        dirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1] / np.sqrt(3),
                         [1, 1, 0] / np.sqrt(2), [0, -1, 1] / np.sqrt(2)])
        site = CoordinationSite(
            "Ni2+", np.zeros(3), dirs[: g.coordination_number] * g.bond_length_mean,
            ("N",) * g.coordination_number,
        )
        v = GeometryVerifier(PDBReference()).verify(_design(site))
        assert v.trust and not v.ood


class TestSelectivity:
    """selectivity_profile re-scores a site as each metal. Geometry discriminates at CN
    extremes but not in the overlapping mid-CN region where most real designs sit."""

    def test_profile_returns_a_verdict_per_metal(self):
        v = GeometryVerifier(PDBReference())
        prof = selectivity_profile(_design(octahedral_site("Ni2+", bond=2.15)), v, ["Ni2+", "Cu2+", "Co2+"])
        assert set(prof) == {"Ni2+", "Cu2+", "Co2+"} and prof["Ni2+"].trust

    def test_cn4_site_does_not_discriminate(self):
        # the honest limit: a CN4 site is in-range for Ni/Cu/Co, so geometry can't select
        v = GeometryVerifier(PDBReference())
        s = octahedral_site("Ni2+", bond=2.15)
        cn4 = CoordinationSite("Ni2+", s.metal_xyz, s.ligand_xyz[:4], s.ligand_elems[:4])
        prof = selectivity_profile(_design(cn4), v, ["Ni2+", "Cu2+", "Co2+"])
        assert all(vd.trust for vd in prof.values())  # trusts for all — geometry gives no selectivity
class TestCSDReference:
    """The CSD/Mogul drop-in: same interface as PDBReference, license-gated data."""

    _CSD_JSON = {
        "Ni2+": {"metal": "Ni2+", "coordination_number": 6, "bond_length_mean": 2.08,
                 "bond_length_std": 0.09, "cn_range": [4, 6], "source": "test fixture"},
    }

    def _ref(self, tmp_path) -> CSDReference:
        p = tmp_path / "csd_reference.json"
        p.write_text(json.dumps(self._CSD_JSON))
        return CSDReference(p)

    def test_plugs_into_verifier_unchanged(self, tmp_path):
        # identical call as PDB/Mock — the verifier never learns which reference it is;
        # an ideal site at the CSD mean (2.08 Å, CN in [4,6]) trusts through it
        v = GeometryVerifier(self._ref(tmp_path)).verify(GOOD)
        assert v.trust and not v.ood

    def test_unknown_metal_raises_with_source(self, tmp_path):
        with pytest.raises(KeyError, match="CSD"):
            self._ref(tmp_path).geometry("Au3+")

    def test_missing_data_raises_actionable_error(self, tmp_path):
        # license-gated: no data file yet ⇒ a clear "build it first" error, not a raw stack trace
        with pytest.raises(FileNotFoundError, match="CSD reference data not found"):
            CSDReference(tmp_path / "nonexistent.json")

    def test_best_reference_falls_back_when_no_metalpdb(self):
        # with no MetalPDB file (the repo default), pick CSD if its data exists, else PDB
        ref = best_reference()
        assert ref.source == ("CSD" if CSD_DATA.exists() else "PDB")
        assert ref.geometry("Ni2+").coordination_number >= 1  # whichever it is, it's usable


class TestMetalPDBReference:
    """The open, licence-free metalloprotein reference — same interface as PDB/CSD."""

    _JSON = {
        "Ni2+": {"metal": "Ni2+", "coordination_number": 5, "bond_length_mean": 2.11,
                 "bond_length_std": 0.12, "cn_range": [4, 6], "source": "test fixture"},
    }

    def _ref(self, tmp_path) -> MetalPDBReference:
        p = tmp_path / "metalpdb_reference.json"
        p.write_text(json.dumps(self._JSON))
        return MetalPDBReference(p)

    def test_plugs_into_verifier_unchanged(self, tmp_path):
        v = GeometryVerifier(self._ref(tmp_path)).verify(GOOD)
        assert isinstance(v.score, float)

    def test_unknown_metal_raises_with_source(self, tmp_path):
        with pytest.raises(KeyError, match="MetalPDB"):
            self._ref(tmp_path).geometry("Au3+")

    def test_missing_data_raises_actionable_error(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="MetalPDB reference data not found"):
            MetalPDBReference(tmp_path / "nonexistent.json")

    def test_best_reference_prefers_metalpdb_when_built(self, tmp_path, monkeypatch):
        # once the open MetalPDB file exists it wins over CSD/PDB — the licence-free default
        p = tmp_path / "metalpdb_reference.json"
        p.write_text(json.dumps(self._JSON))
        monkeypatch.setattr(_ref_mod, "METALPDB_DATA", p)
        ref = best_reference()
        assert ref.source == "MetalPDB" and ref.geometry("Ni2+").bond_length_mean == 2.11


class TestMetalHawk:
    """The open ANN geometry-distortion oracle — a learned CN/geometry classifier, pluggable
    via a scorer callback (the heavy inference lives in scripts/metalhawk_score.py)."""

    def _design(self, n: int) -> BinderDesign:
        s = octahedral_site("Ni2+")
        site = CoordinationSite("Ni2+", s.metal_xyz, s.ligand_xyz[:n], s.ligand_elems[:n])
        return BinderDesign("SEQ", site, generator="t", generator_confidence=0.5, source="x.pdb")

    def test_confident_cn_match_is_trusted(self):
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(6, "octahedral", 0.9)).verify(self._design(6))
        assert v.trust and not v.ood and v.score > 0.8

    def test_small_cn_mismatch_is_weak(self):
        # Δcn == 1: MetalHawk confidently sees a mildly different geometry ⇒ a distortion signal, weak
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(5, "square pyramidal", 0.9)).verify(self._design(6))
        assert not v.trust and not v.ood and v.score < 0.6

    def test_confident_large_mismatch_defers_off_manifold(self):
        # Δcn >= ood_cn_gap at high confidence: MetalHawk grossly contradicts the physical shell
        # ⇒ confidently off its training manifold, so it abstains rather than emit a spurious weak.
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(4, "tetrahedral", 0.99)).verify(self._design(6))
        assert v.ood and not v.trust and "off-manifold" in v.reason

    def test_low_confidence_defers(self):
        # the ANN itself is unsure ⇒ off its training manifold, not judgeable
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(6, "octahedral", 0.2)).verify(self._design(6))
        assert v.ood and not v.trust

    @pytest.mark.parametrize("scorer", [lambda _d: None, _boom], ids=["no-prediction", "scorer-raises"])
    def test_defers_without_a_prediction(self, scorer):
        v = MetalHawkVerifier(scorer).verify(self._design(6))
        assert v.ood and not v.trust and "MetalHawk" in v.reason

    def test_score_provider_keys_by_source(self):
        pred = MetalHawkPrediction(6, "octahedral", 0.9)
        provide = score_provider({"x.pdb": pred})
        assert provide(self._design(6)) is pred  # source "x.pdb" → its prediction
        other = self._design(6)
        other.source = "unlisted.pdb"
        assert provide(other) is None  # no prediction for an unlisted structure

    def test_load_predictions_round_trips(self, tmp_path):
        # the loader that closes the metalhawk_score.py JSON → verifier loop
        p = tmp_path / "scores.json"
        p.write_text(json.dumps({"a.pdb": {"coordination_number": 4, "geometry": "tetrahedral", "confidence": 0.8}}))
        assert load_predictions(p) == {"a.pdb": MetalHawkPrediction(4, "tetrahedral", 0.8)}


class TestEmptySite:
    def test_no_donors_defers_with_finite_score(self, verifier):
        """A site with no coordinating atoms (CN=0) is the worst case, not NaN."""
        empty = _design(CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ()))
        v = verifier.verify(empty)
        assert v.score == 0.0 and not v.trust and v.ood
        assert rank([empty, GOOD], verifier)[0][0] is GOOD  # NaN would break this sort
