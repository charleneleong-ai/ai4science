import importlib.util
import json
import re
from pathlib import Path

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
from touchstone.geometry import reference as ref_mod
from touchstone.geometry.metalhawk import load_predictions, score_provider
from touchstone.geometry.parse import SOLVENT_RESIDUES, cif_atoms, coordination_site, pdb_atoms
from touchstone.geometry.reference import (
    CSD_DATA,
    METALPDB_DATA,
    ReferenceDistribution,
    best_reference,
)

# the reference builder is a script (not importable as a package) — load it by path to unit-test
# the domain guards (solvent exclusion, sparse-metal floor) without hitting the MetalPDB API
spec = importlib.util.spec_from_file_location(
    "build_metalpdb_reference", Path(__file__).parent.parent / "scripts" / "build_metalpdb_reference.py"
)
bmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bmp)

FIX = Path(__file__).parent / "fixtures"


def boom(_design):  # a scorer that fails — module-level so it's not a lambda
    raise RuntimeError("metalhawk exploded")


def design(site: CoordinationSite, conf: float = 0.7) -> BinderDesign:
    return BinderDesign("SEQ", site, generator="test", generator_confidence=conf)


def site_at(metal: str, donors: list[float]) -> CoordinationSite:
    """A site whose donors sit at the given distances (Å), spread over a golden-spiral sphere.
    The geometry tier scores bond lengths + CN, so the directions only need to be distinct."""
    i = np.arange(len(donors)) + 0.5
    phi = np.arccos(1 - 2 * i / len(donors))
    theta = np.pi * (1 + 5**0.5) * i
    dirs = np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
    return CoordinationSite(
        metal=metal,
        metal_xyz=np.zeros(3),
        ligand_xyz=dirs * np.array(donors)[:, None],
        ligand_elems=("N",) * len(donors),
    )


@pytest.fixture(scope="module")
def real_sites() -> dict[str, list[list[float]]]:
    """Real metalloprotein first-shell donor distances (MetalPDB), by metal — the ground truth a
    geometry prior must not reject."""
    raw = json.loads((FIX / "metalpdb_real_sites.json").read_text())
    return {metal: [s["donors"] for s in sites] for metal, sites in raw.items()}


@pytest.fixture
def verifier() -> GeometryVerifier:
    return GeometryVerifier(MockReference())


# A near-ideal Ni2+ site vs. two ways to be wrong: stretched bonds, and wrong coordination number.
GOOD = design(octahedral_site("Ni2+", bond=2.09))
STRAINED = design(octahedral_site("Ni2+", bond=2.6))
WRONG_CN = design(
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
        low = verifier.verify(design(site, conf=0.01))
        high = verifier.verify(design(site, conf=0.99))
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
            verifier.verify(design(octahedral_site("Au3+")))


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
        v = GeometryVerifier(PDBReference()).verify(design(site))
        assert v.trust and not v.ood


class TestSelectivity:
    """selectivity_profile re-scores a site as each metal. Geometry discriminates at CN
    extremes but not in the overlapping mid-CN region where most real designs sit."""

    def test_profile_returns_a_verdict_per_metal(self):
        v = GeometryVerifier(PDBReference())
        prof = selectivity_profile(design(octahedral_site("Ni2+", bond=2.15)), v, ["Ni2+", "Cu2+", "Co2+"])
        assert set(prof) == {"Ni2+", "Cu2+", "Co2+"} and prof["Ni2+"].trust

    def test_cn4_site_does_not_discriminate(self):
        # the honest limit: a CN4 site is in-range for Ni/Cu/Co, so geometry can't select
        v = GeometryVerifier(PDBReference())
        s = octahedral_site("Ni2+", bond=2.15)
        cn4 = CoordinationSite("Ni2+", s.metal_xyz, s.ligand_xyz[:4], s.ligand_elems[:4])
        prof = selectivity_profile(design(cn4), v, ["Ni2+", "Cu2+", "Co2+"])
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

    def test_best_reference_falls_back_when_no_metalpdb(self, tmp_path, monkeypatch):
        # MetalPDB now ships, so force it absent to exercise the fallback: CSD if built, else PDB
        monkeypatch.setattr(ref_mod, "METALPDB_DATA", tmp_path / "absent.json")
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
        monkeypatch.setattr(ref_mod, "METALPDB_DATA", p)
        ref = best_reference()
        assert ref.source == "MetalPDB" and ref.geometry("Ni2+").bond_length_mean == 2.11


def trust_rate(shells: list[list[float]], metal: str, ref: ReferenceDistribution) -> float:
    """Fraction of real donor shells the geometry tier trusts under `ref`."""
    v = GeometryVerifier(ref)
    verdicts = [v.verify(design(site_at(metal, shell))) for shell in shells]
    return sum(x.trust for x in verdicts) / len(verdicts)


class TestShippedPriorCalibration:
    """The prior that ships, checked against real, experimentally-determined metalloprotein sites
    (tests/fixtures/metalpdb_real_sites.json — 360 MetalPDB sites, protein first-shell donors,
    regenerated by scripts/build_real_site_fixture.py through the builder's own filters).

    This is in-sample — the fixture is a subsample of the prior's own source — so it is a *coverage*
    check, not held-out validation: does the prior cover the domain it claims to describe? That is a
    low bar, and the point is that the CSD prior fails it (trusting 41% of real Ni sites, deferring
    19.5% as off-manifold) while the shipped one clears it.
    See docs/experiments/2026-07-13-geometry-prior-wrong-domain.md."""

    @pytest.mark.parametrize("metal", ["Ni2+", "Cu2+", "Co2+"])
    def test_shipped_prior_trusts_real_metalloprotein_sites(self, metal, real_sites):
        # a 2σ gate on a correctly-centred prior should admit ~95% of real sites (Gaussian), and
        # does: Ni 99%, Cu 98%, Co 96%. This is the calibration claim and the guard on the tier.
        rate = trust_rate(real_sites[metal], metal, best_reference())
        assert rate >= 0.90, f"shipped prior trusts only {rate:.0%} of real {metal} protein sites"

    @pytest.mark.parametrize(
        "shell,why",
        [([3.5, 3.5, 3.5, 3.5], "donors 3.5 Å out"), ([2.19], "a single donor"), ([2.19] * 12, "CN 12")],
        ids=["stretched", "cn1", "cn12"],
    )
    def test_implausible_sites_are_still_rejected(self, shell, why):
        # the other half of the calibration claim. Sensitivity without specificity is free: a prior
        # asserting Ni-donor bonds are 9.9 Å long would trust every real site too. Widening the prior
        # must not become a way to pass the tests above.
        v = GeometryVerifier(best_reference()).verify(design(site_at("Ni2+", shell)))
        assert not v.trust, f"shipped prior trusts a Ni site with {why}"

    def test_octahedral_copper_is_rejected(self):
        """Protein Cu2+ is 3-5 coordinate (type-1/type-2 copper is not octahedral); octahedral Cu is
        small-molecule/aqua geometry. The CSD prior rates this site `plausible (0.0σ)` — a *perfect*
        score — because its cn_range runs to 6. That is how a wrong-domain prior hides from a
        real-site coverage check: it fails on specificity, not sensitivity. It also means a Cu RLVR
        geometry reward on the CSD prior would have rewarded octahedral copper."""
        cu6 = design(site_at("Cu2+", [2.12] * 6))
        assert GeometryVerifier(CSDReference()).verify(cu6).trust  # the bug, pinned
        assert not GeometryVerifier(best_reference()).verify(cu6).trust

    def test_the_csd_prior_would_fail_this(self, real_sites):
        # the regression, stated as a test: the previous default rejects real nickel biochemistry.
        # If this ever starts passing, the CSD data changed domain and the finding needs revisiting.
        rate = trust_rate(real_sites["Ni2+"], "Ni2+", CSDReference())
        assert rate < 0.60, f"CSD unexpectedly trusts {rate:.0%} of real Ni sites — re-check the prior audit"

    def test_default_prior_is_the_protein_domain_one(self):
        assert best_reference().source == "MetalPDB"

    def test_every_shipped_prior_clears_the_site_floor(self):
        # ties the artifact back to the constant — deleting MIN_SITES from the builder must not leave
        # a thin prior in the bundle. The non-empty check matters too: best_reference() selects on
        # file existence, so an empty table would be "present" and KeyError on every metal.
        table = json.loads(METALPDB_DATA.read_text())
        assert {"Ni2+", "Cu2+", "Co2+"} <= set(table)
        for metal, g in table.items():
            n = int(re.search(r"(\d+) sites", g["source"]).group(1))
            assert n >= bmp.MIN_SITES, f"{metal} prior built from only {n} sites"

    def test_builder_refuses_to_write_an_empty_table(self, monkeypatch):
        # an offline run must not clobber the shipped prior with {} — "present but empty" is worse
        # than absent, because the preference order selects on existence
        monkeypatch.setattr(bmp, "fetch", lambda _symbol: [])
        monkeypatch.setattr(bmp, "OUT", Path("/nonexistent/must-not-be-written.json"))
        with pytest.raises(SystemExit):
            bmp.main()


class TestSolventExclusion:
    """The parser must count the same donor set the prior was measured on: protein only. A design
    carries no waters, but a crystal structure or a co-fold does — and counting its water oxygens
    against a solvent-free prior is the domain error of
    docs/experiments/2026-07-13-geometry-prior-wrong-domain.md, mirrored."""

    FIXTURES = [FIX / "nickel_site_with_waters.pdb", FIX / "nickel_site_with_waters.cif"]

    @pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.suffix)
    def test_water_donors_are_not_counted(self, path):
        site = coordination_site(path, "NI", "Ni2+", 2.8)
        assert site.coordination_number == 3
        assert sorted(site.ligand_elems) == ["N", "O", "S"]  # His/Asp/Cys, not the two waters

    @pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.suffix)
    def test_the_waters_are_inside_the_cutoff(self, path):
        # without this the test above could pass for the wrong reason (waters simply too far away).
        # Both fixture waters sit ~2.1 Å from the metal — an element-only parser would see CN=5.
        waters = [
            xyz for el, res, xyz in (pdb_atoms if path.suffix == ".pdb" else cif_atoms)(path.read_text())
            if res in SOLVENT_RESIDUES and el == "O"
        ]
        assert len(waters) == 2
        assert all(np.linalg.norm(w) <= 2.8 for w in waters)


class TestMetalPDBReferenceBuilder:
    """Solvent exclusion, on synthetic records (no network). The site floor is checked against
    the shipped artifact in TestShippedPriorCalibration."""

    def record(self, *residues: str) -> list[dict]:
        # one Ni site whose donors come from the given residues, all at 2.20 Å
        return [{"metals": [{"symbol": "Ni", "ligands": [
            {"residue": r, "donors": [{"symbol": "O", "distance": 2.20}]} for r in residues
        ]}]}]

    def test_solvent_donors_are_excluded(self):
        # 3 protein donors + 2 waters: a design carries no solvent, so the prior must count 3.
        # Counting the waters is what pushes Ni/Co/Mn modal CN 6→4 and marks a good site with an
        # open, water-fillable face as under-coordinated.
        rec = self.record("HIS", "ASP", "CYS", "HOH", "HOH")
        assert [len(s) for s in bmp.donor_shells(rec, "Ni")] == [3]


class TestMetalHawk:
    """The open ANN geometry-distortion oracle — a learned CN/geometry classifier, pluggable
    via a scorer callback (the heavy inference lives in scripts/metalhawk_score.py)."""

    def design(self, n: int) -> BinderDesign:
        s = octahedral_site("Ni2+")
        site = CoordinationSite("Ni2+", s.metal_xyz, s.ligand_xyz[:n], s.ligand_elems[:n])
        return BinderDesign("SEQ", site, generator="t", generator_confidence=0.5, source="x.pdb")

    def test_confident_cn_match_is_trusted(self):
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(6, "octahedral", 0.9)).verify(self.design(6))
        assert v.trust and not v.ood and v.score > 0.8

    def test_small_cn_mismatch_is_weak(self):
        # Δcn == 1: MetalHawk confidently sees a mildly different geometry ⇒ a distortion signal, weak
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(5, "square pyramidal", 0.9)).verify(self.design(6))
        assert not v.trust and not v.ood and v.score < 0.6

    def test_confident_large_mismatch_defers_off_manifold(self):
        # Δcn >= ood_cn_gap at high confidence: MetalHawk grossly contradicts the physical shell
        # ⇒ confidently off its training manifold, so it abstains rather than emit a spurious weak.
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(4, "tetrahedral", 0.99)).verify(self.design(6))
        assert v.ood and not v.trust and "off-manifold" in v.reason

    def test_low_confidence_defers(self):
        # the ANN itself is unsure ⇒ off its training manifold, not judgeable
        v = MetalHawkVerifier(lambda _d: MetalHawkPrediction(6, "octahedral", 0.2)).verify(self.design(6))
        assert v.ood and not v.trust

    @pytest.mark.parametrize("scorer", [lambda _d: None, boom], ids=["no-prediction", "scorer-raises"])
    def test_defers_without_a_prediction(self, scorer):
        v = MetalHawkVerifier(scorer).verify(self.design(6))
        assert v.ood and not v.trust and "MetalHawk" in v.reason

    def test_score_provider_keys_by_source(self):
        pred = MetalHawkPrediction(6, "octahedral", 0.9)
        provide = score_provider({"x.pdb": pred})
        assert provide(self.design(6)) is pred  # source "x.pdb" → its prediction
        other = self.design(6)
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
        empty = design(CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ()))
        v = verifier.verify(empty)
        assert v.score == 0.0 and not v.trust and v.ood
        assert rank([empty, GOOD], verifier)[0][0] is GOOD  # NaN would break this sort
