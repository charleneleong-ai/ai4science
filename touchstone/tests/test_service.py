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
        assert all("label" in v and "score" in v for v in r["verifiers"].values())
        assert {"mogul", "cofold", "expression", "thermostability"} <= set(r["not_run"])  # full stack advertised

    def test_cif_goes_through_the_same_path(self):
        r = verify_structure(FIX / "boltzgen_nickel_design.cif", "Ni2+")
        assert r["metal"] == "Ni2+" and r["consensus"] in {"trust", "weak", "defer"}

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
