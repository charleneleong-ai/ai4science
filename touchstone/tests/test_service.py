import json
from pathlib import Path

from typer.testing import CliRunner

from touchstone.cli import app
from touchstone.service import verify_structure

FIX = Path(__file__).parent / "fixtures"
PACKED = FIX / "ligmpnn_nickel_packed.pdb"
runner = CliRunner()


class TestVerifyStructure:
    def test_returns_per_verifier_and_consensus(self):
        r = verify_structure(PACKED, "Ni2+")
        assert set(r["verifiers"]) == {"geometry", "bond_valence"}  # lightweight default, no GPU
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
            "geometry", "bond_valence", "mogul", "mlip", "mlip_md", "cofold", "expression", "thermostability"
        ]
        assert by["geometry"]["status"] == "ran" and "strain_sigma" in by["geometry"]["metrics"]
        assert by["bond_valence"]["status"] == "ran" and "bvs" in by["bond_valence"]["metrics"]
        assert by["mogul"]["status"] == "needs_input"  # CSD licence
        assert by["mlip"]["status"] == "needs_input"  # needs deep=True + a GPU

    def test_stress_adds_a_robustness_map_only_when_requested(self):
        assert "stress" not in verify_structure(PACKED, "Ni2+")
        r = verify_structure(PACKED, "Ni2+", stress=True)
        assert set(r["stress"]) == {"neutral", "leachate", "low_pH"}  # the operating-envelope conditions
        assert all(v["label"] in {"trust", "weak", "defer"} for v in r["stress"].values())

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
