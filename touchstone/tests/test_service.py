import json
from pathlib import Path

from typer.testing import CliRunner

from touchstone import ExpressionSignals, MetalHawkPrediction, ThermostabilitySignal
from touchstone.cli import app
from touchstone.cofold import cif_provider
from touchstone.service import verify_structure

FIX = Path(__file__).parent / "fixtures"
PACKED = FIX / "ligmpnn_nickel_packed.pdb"
runner = CliRunner()


class TestVerifyStructure:
    def test_returns_per_verifier_and_consensus(self):
        r = verify_structure(PACKED, "Ni2+")
        # the always-on pure-python tiers (no GPU): lengths, valence, symmetry, shape
        assert set(r["verifiers"]) == {"geometry", "bond_valence", "coord_symmetry", "coord_geometry", "precedent"}
        assert r["consensus"] in {"trust", "weak", "defer"}
        assert r["coordination_number"] >= 1
        assert r["reference"] in {"CSD", "PDB"}  # which geometry prior backed the z-score
        assert all("label" in v and "score" in v for v in r["verifiers"].values())
        assert {"mogul", "cofold", "expression", "thermostability"} <= set(r["not_run"])  # full stack advertised

    def test_cif_goes_through_the_same_path(self):
        r = verify_structure(FIX / "boltzgen_nickel_design.cif", "Ni2+")
        assert r["metal"] == "Ni2+" and r["consensus"] in {"trust", "weak", "defer"}

    def test_stack_lists_every_tier_with_status_and_metrics(self):
        r = verify_structure(PACKED, "Ni2+")
        by = {s["stage"]: s for s in r["stack"]}
        # the full stack appears in cost order, even tiers that didn't run on a bare structure
        assert [s["stage"] for s in r["stack"]] == [
            "geometry", "bond_valence", "coord_symmetry", "coord_geometry", "precedent", "metalhawk",
            "mogul", "mlip", "mlip_md", "trs", "selectivity", "cofold", "expression", "thermostability",
        ]
        assert by["geometry"]["status"] == "ran" and "strain_sigma" in by["geometry"]["metrics"]
        assert by["bond_valence"]["status"] == "ran" and "bvs" in by["bond_valence"]["metrics"]
        assert by["coord_symmetry"]["status"] == "ran" and "nvecsum" in by["coord_symmetry"]["metrics"]
        assert by["coord_geometry"]["status"] == "ran" and "angle_rmsd_deg" in by["coord_geometry"]["metrics"]
        assert by["metalhawk"]["status"] == "needs_input"  # MetalHawk predictions (open, no licence)
        assert by["mogul"]["status"] == "needs_input"  # CSD licence
        assert by["mlip"]["status"] == "needs_input"  # needs deep=True + a GPU
        assert by["selectivity"]["status"] == "needs_input"  # needs deep=True + selectivity_metals

    def test_stress_adds_a_robustness_map_only_when_requested(self):
        assert "stress" not in verify_structure(PACKED, "Ni2+")
        r = verify_structure(PACKED, "Ni2+", stress=True)
        assert set(r["stress"]) == {"neutral", "leachate", "low_pH"}  # the operating-envelope conditions
        assert all(v["label"] in {"trust", "weak", "defer"} for v in r["stress"].values())

    def test_cofold_provider_adds_the_crosscheck_tier(self):
        # provider predicts the design's own structure ⇒ the independent site agrees ⇒ trust.
        provider = cif_provider({str(PACKED): str(PACKED)}, metal_atom="NI", metal="Ni2+")
        r = verify_structure(PACKED, "Ni2+", cofold_provider=provider)
        by = {s["stage"]: s for s in r["stack"]}
        assert by["cofold"]["status"] == "ran"  # was needs_input without a provider
        assert r["verifiers"]["cofold"]["label"] == "trust"

    def test_metalhawk_scorer_adds_the_tier(self):
        # scorer echoes the design's CN ⇒ MetalHawk agrees ⇒ the tier runs and trusts.
        scorer = lambda d: MetalHawkPrediction(d.site.coordination_number, "octahedral", 0.9)
        r = verify_structure(PACKED, "Ni2+", metalhawk_scorer=scorer)
        by = {s["stage"]: s for s in r["stack"]}
        assert by["metalhawk"]["status"] == "ran"  # was needs_input without a scorer
        assert r["verifiers"]["metalhawk"]["label"] == "trust"

    def test_opt_in_kwargs_wire_their_tiers(self):
        # the sequence-keyed seams flip their tiers from needs_input to ran (precedent runs by default)
        r = verify_structure(
            PACKED, "Ni2+",
            expression_scorer=lambda d: ExpressionSignals(8.0, 0.7),
            thermostability_predictor=lambda d: ThermostabilitySignal(65.0),
        )
        by = {s["stage"]: s["status"] for s in r["stack"]}
        assert by["expression"] == "ran" and by["thermostability"] == "ran"

    def test_selectivity_metals_wires_it_into_the_deep_family(self):
        # off without a metals panel; with it (+ deep) it joins mlip/mlip_md. calc=None forces the
        # no-backend path, so it's "skipped" (proving it's wired) rather than "needs_input".
        default = {s["stage"]: s["status"] for s in verify_structure(PACKED, "Ni2+")["stack"]}
        assert default["selectivity"] == "needs_input"
        wired = verify_structure(PACKED, "Ni2+", deep=True, calc=None, selectivity_metals=("Ni2+", "Cu2+", "Co2+"))
        assert {s["stage"]: s["status"] for s in wired["stack"]}["selectivity"] == "skipped"

    def test_deep_without_backend_degrades_gracefully(self):
        # no GPU/mace here ⇒ mlip + mlip_md are skipped, consensus still decided by geometry+BV
        r = verify_structure(PACKED, "Ni2+", deep=True)
        assert "skipped" in r["verifiers"]["mlip"] and "skipped" in r["verifiers"]["mlip_md"]
        assert r["consensus"] in {"trust", "weak", "defer"}


class TestCLI:
    def test_human_output_shows_consensus(self):
        res = runner.invoke(app, ["verify", str(PACKED), "--metal", "Ni2+"])
        assert res.exit_code == 0 and "consensus" in res.stdout.lower()

    def test_json_output_is_parseable(self):
        res = runner.invoke(app, ["verify", str(PACKED), "--json"])
        assert res.exit_code == 0
        data = json.loads(res.stdout)
        assert data["metal"] == "Ni2+" and "consensus" in data and "verifiers" in data
