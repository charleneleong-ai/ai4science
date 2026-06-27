from pathlib import Path

from typer.testing import CliRunner

from touchstone.cli import app
from touchstone.reward import best_of_n, rank_structures, reward_from_result

FIX = Path(__file__).parent / "fixtures"
DESIGNS = [FIX / "ligmpnn_nickel_packed.pdb", FIX / "rfaa_nickel_sample_0.pdb", FIX / "boltzgen_nickel_design.cif"]
runner = CliRunner()


def _result(consensus: str, *scores: float) -> dict:
    return {"consensus": consensus, "verifiers": {f"v{i}": {"score": s} for i, s in enumerate(scores)}}


class TestRewardFromResult:
    def test_consensus_scales_the_reward(self):
        # same scores, trust > weak > defer
        trust = reward_from_result(_result("trust", 0.8, 0.6))
        weak = reward_from_result(_result("weak", 0.8, 0.6))
        defer = reward_from_result(_result("defer", 0.8, 0.6))
        assert trust > weak > defer
        assert defer == 0.0 and trust == 0.7  # mean(0.8,0.6) × 1.0

    def test_no_scores_is_zero(self):
        assert reward_from_result({"consensus": "defer", "verifiers": {}}) == 0.0


class TestRankStructures:
    def test_sorted_best_first_with_rewards(self):
        ranked = rank_structures(DESIGNS, "Ni2+")
        rewards = [r["reward"] for r in ranked]
        assert rewards == sorted(rewards, reverse=True)
        assert all(0.0 <= r <= 1.0 for r in rewards)

    def test_unparseable_structure_scores_zero(self, tmp_path):
        junk = tmp_path / "no_metal.pdb"
        junk.write_text("ATOM      1  N   HIS A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
        ranked = rank_structures([junk], "Ni2+")
        assert ranked[0]["reward"] == 0.0 and "error" in ranked[0]

    def test_best_of_n_takes_the_top(self):
        top = best_of_n(DESIGNS, "Ni2+", n=1)
        assert len(top) == 1 and top[0]["reward"] == max(r["reward"] for r in rank_structures(DESIGNS, "Ni2+"))


def test_cli_rank_runs():
    res = runner.invoke(app, ["rank", *[str(d) for d in DESIGNS], "--metal", "Ni2+"])
    assert res.exit_code == 0 and "rank" in res.stdout.lower()
