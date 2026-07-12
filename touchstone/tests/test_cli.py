import json
from pathlib import Path

from typer.testing import CliRunner

from touchstone.cli import app

runner = CliRunner()
FIXTURE = str(Path(__file__).parent / "fixtures" / "ideal_nickel_site.pdb")


def _verify(*args) -> dict:
    result = runner.invoke(app, ["verify", FIXTURE, "--metal", "Ni2+", "--json", *args])
    assert result.exit_code == 0, result.output
    return json.loads(result.output)


def _status(r: dict, stage: str) -> str:
    return next(s["status"] for s in r["stack"] if s["stage"] == stage)


class TestVerifyCLI:
    """Tier flags on `touchstone verify` — precedent on by default, the others opt-in (site/path-based)."""

    def test_precedent_on_by_default(self):
        r = _verify()  # open MetalPDB precedent runs out of the box
        assert _status(r, "precedent") == "ran" and "precedent(s)" in r["verifiers"]["precedent"]["reason"]

    def test_no_precedent_disables_it(self):
        r = _verify("--no-precedent")
        assert "precedent" not in {s["stage"] for s in r["stack"]}

    def test_metalhawk_scores_enables_the_tier(self, tmp_path):
        scores = tmp_path / "mh.json"
        scores.write_text(json.dumps({FIXTURE: {"coordination_number": 6, "geometry": "octahedral", "confidence": 0.9}}))
        assert _status(_verify("--metalhawk-scores", str(scores)), "metalhawk") == "ran"

    def test_sequence_enables_the_sequence_keyed_tiers(self, tmp_path):
        # expression / thermostability key by sequence — --sequence bridges the score files to the design
        seq = "MKVLAAA"
        expr = tmp_path / "expr.json"
        expr.write_text(json.dumps({seq: {"pseudo_perplexity": 8.0, "solubility": 0.7}}))
        tm = tmp_path / "tm.json"
        tm.write_text(json.dumps({seq: 65.0}))
        r = _verify("--sequence", seq, "--expression-scores", str(expr), "--thermostability-scores", str(tm))
        assert (_status(r, "expression"), r["verifiers"]["expression"]["label"]) == ("ran", "trust")
        assert (_status(r, "thermostability"), r["verifiers"]["thermostability"]["label"]) == ("ran", "trust")


class TestRankCLI:
    def test_rank_precedent_flag_folds_into_the_reward(self):
        result = runner.invoke(app, ["rank", FIXTURE, "--metal", "Ni2+", "--precedent", "--json"])
        assert result.exit_code == 0, result.output
        ranked = json.loads(result.output)
        assert ranked and "reward" in ranked[0]
