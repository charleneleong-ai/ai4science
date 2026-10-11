"""What gets logged to W&B: the funnel's own consistency, and the assets the logger collects.

`log()` itself needs wandb and a network, so it is not exercised here. Everything it reads --
the funnel, the ranked designs, the structure list, the scan CSVs -- is pure, and that is what
decides whether the run is honest.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from conftest import MODEL_0

from tnf_binder.tracking import (FUNNEL, FUNNELS, LEVERS, SUBMISSION2_FUNNEL,
                                 CampaignLogger, funnel_surviving)


def build_assets(root: Path, designs: dict[str, dict[str, object]]) -> Path:
    """An assets/ tree shaped like the campaign's: boltz_camp/<design>/{confidence,model}."""
    for name, payload in designs.items():
        directory = root / "boltz_camp" / name
        directory.mkdir(parents=True)
        (directory / f"confidence_{name}_model_0.json").write_text(json.dumps(payload))
        (directory / f"{name}_model_0.cif").write_text("data_synthetic\n")
    return root


def weakened(payload: dict[str, object], value: float) -> dict[str, object]:
    """The same sample with the binder's pair ipTM flattened to `value`."""
    pairs = {k: dict(v) for k, v in payload["pair_chains_iptm"].items()}
    for target in ("1", "2", "3"):
        pairs["0"][target] = value
        pairs[target]["0"] = value
    return {**payload, "pair_chains_iptm": pairs}


# all three share MODEL_0's global ipTM (0.771) and differ only in the binder's pair ipTM,
# so the ranking under test is the interface and nothing else.
STRONG, MIDDLING, WEAK = 0.62, 0.40, 0.05


@pytest.fixture
def logger(tmp_path) -> CampaignLogger:
    assets = build_assets(tmp_path, {
        "strong": weakened(MODEL_0, STRONG),
        "weak": weakened(MODEL_0, WEAK),
        "middling": weakened(MODEL_0, MIDDLING),
    })
    return CampaignLogger(assets)


class TestFunnel:
    """The funnel is hand-maintained, so its arithmetic is asserted rather than assumed."""

    def test_each_stage_enters_with_the_previous_stage_survivors(self):
        for (_, entering, _), (previous, _, surviving) in zip(FUNNEL[1:], FUNNEL):
            assert entering == surviving, f"{previous!r} -> entering {entering}"

    def test_no_stage_passes_more_than_it_received(self):
        assert all(surviving <= entering for _, entering, surviving in FUNNEL)

    @pytest.mark.parametrize("stage,expected", [("submitted", 3), ("novelty", 2),
                                                ("generated", 57)])
    def test_lookup_by_fragment(self, stage, expected):
        assert funnel_surviving(stage) == expected

    def test_unknown_stage_raises(self):
        with pytest.raises(KeyError):
            funnel_surviving("expression")


class TestSubmission2Funnel:
    """The second submission has its own denominators and no novelty verdict yet."""

    def test_arithmetic_holds_for_every_registered_funnel(self):
        for name, funnel in FUNNELS.items():
            assert all(v <= e for _, e, v in funnel), name

    def test_summary_omits_a_novelty_count_it_cannot_have(self, tmp_path):
        build_assets(tmp_path, {"a": MODEL_0})
        summary = CampaignLogger(tmp_path, funnel=SUBMISSION2_FUNNEL).summary()
        assert summary["designs_submitted"] == 12
        assert "designs_passing_novelty" not in summary

    def test_campaign_summary_still_reports_novelty(self, tmp_path):
        build_assets(tmp_path, {"a": MODEL_0})
        assert CampaignLogger(tmp_path).summary()["designs_passing_novelty"] == 2


class TestLevers:
    """Only one lever moved the interface; the table must keep saying so."""

    def test_exactly_one_lever_had_an_effect(self):
        moved = [name for name, _, delta in LEVERS if abs(delta) >= 0.05]
        assert moved == ["target stoichiometry"]


class TestCampaignLogger:
    """What the logger collects before wandb is involved."""

    def test_designs_are_ranked_best_interface_first(self, logger):
        assert [d.name for d in logger.designs] == ["strong", "middling", "weak"]

    def test_structures_follow_the_design_ranking(self, logger):
        assert [p.parent.name for p in logger.structures] == ["strong", "middling", "weak"]

    def test_structures_are_capped(self, tmp_path):
        assets = build_assets(tmp_path, {f"d{i}": MODEL_0 for i in range(5)})
        assert len(CampaignLogger(assets, max_molecules=2).structures) == 2

    def test_summary_counts_come_from_the_funnel(self, logger):
        summary = logger.summary()
        assert summary["designs_folded"] == 3
        assert summary["designs_submitted"] == funnel_surviving("submitted")
        assert summary["designs_passing_novelty"] == funnel_surviving("novelty")

    def test_summary_reports_the_best_design_and_the_inflation(self, logger):
        summary = logger.summary()
        assert summary["best_binder_to_target_iptm"] == pytest.approx(STRONG, abs=1e-3)
        # global ipTM (0.771) minus each interface, averaged -- the artefact, quantified
        expected = 0.771 - (STRONG + MIDDLING + WEAK) / 3
        assert summary["mean_global_iptm_inflation"] == pytest.approx(expected, abs=2e-3)

    def test_read_csv_splits_header_from_rows(self, tmp_path):
        path = tmp_path / "scan.csv"
        path.write_text("parent,pka,verdict\nhelix_0,6.66,ideal\nhelix_0,3.22,dead\n")
        columns, rows = CampaignLogger.read_csv(path)
        assert columns == ["parent", "pka", "verdict"]
        assert rows == [["helix_0", "6.66", "ideal"], ["helix_0", "3.22", "dead"]]

    def test_logs_every_csv_beside_the_predictions(self, logger):
        """Naming individual scans meant a newly-added scan was silently never logged."""
        (logger.assets / "pka_deep.csv").write_text("parent,pka\nhelix_0,6.66\n")
        (logger.assets / "evidence.csv").write_text("name,xr\ntnf2_01,0.79\n")
        columns, rows = CampaignLogger.read_csv(logger.assets / "pka_deep.csv")
        assert columns == ["parent", "pka"] and rows == [["helix_0", "6.66"]]
        assert {p.stem for p in logger.assets.glob("*.csv")} == {"pka_deep", "evidence"}

    def test_predictions_directory_is_selectable(self, tmp_path):
        """Submission 2 folded into its own directory beside the campaign's."""
        build_assets(tmp_path, {"a": MODEL_0})
        other = tmp_path / "boltz_sub2" / "b"
        other.mkdir(parents=True)
        (other / "confidence_b_model_0.json").write_text(json.dumps(weakened(MODEL_0, 0.7)))
        assert [d.name for d in CampaignLogger(tmp_path).designs] == ["a"]
        assert [d.name for d in
                CampaignLogger(tmp_path, predictions_dir="boltz_sub2").designs] == ["b"]

    def test_read_csv_of_an_empty_file(self, tmp_path):
        path = tmp_path / "empty.csv"
        path.write_text("")
        assert CampaignLogger.read_csv(path) == ([], [])

    def test_an_empty_campaign_yields_nothing_rather_than_raising(self, tmp_path):
        logger = CampaignLogger(tmp_path)
        assert logger.designs == []
        assert logger.structures == []
        assert logger.summary()["designs_folded"] == 0
