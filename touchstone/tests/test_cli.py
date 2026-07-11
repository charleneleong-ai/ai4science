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
    """The opt-in tier flags on `touchstone verify` (site/path-based tiers only)."""

    def test_precedent_flag_enables_the_open_tier(self):
        r = _verify("--precedent")
        assert _status(r, "precedent") == "ran" and "precedent(s)" in r["verifiers"]["precedent"]["reason"]

    def test_precedent_off_by_default(self):
        assert _status(_verify(), "precedent") == "needs_input"

    def test_metalhawk_scores_enables_the_tier(self, tmp_path):
        scores = tmp_path / "mh.json"
        scores.write_text(json.dumps({FIXTURE: {"coordination_number": 6, "geometry": "octahedral", "confidence": 0.9}}))
        assert _status(_verify("--metalhawk-scores", str(scores)), "metalhawk") == "ran"


class TestRankCLI:
    def test_rank_precedent_flag_folds_into_the_reward(self):
        result = runner.invoke(app, ["rank", FIXTURE, "--metal", "Ni2+", "--precedent", "--json"])
        assert result.exit_code == 0, result.output
        ranked = json.loads(result.output)
        assert ranked and "reward" in ranked[0]
