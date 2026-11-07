"""The two YAML files the report and the artifact page both render from.

`verdicts.yaml` shipped for a week unparseable -- a bare `: ` inside a list item makes YAML read
the line as a mapping -- because it had only ever been checked with a regex. Anything that two
deliverables render from gets loaded by a parser in CI.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import build_artifact
import build_wandb_report
import common

DOCS = Path(__file__).resolve().parent.parent / "docs" / "experiments" / "pmhc"
VERDICTS_PATH = DOCS / "verdicts.yaml"
CAMPAIGN_PATH = DOCS / "campaign.yaml"
VERDICT_VALUES = set(build_artifact.VERDICT_ORDER)  # what the page can actually colour


@pytest.fixture(scope="module")
def verdicts() -> dict[str, object]:
    return yaml.safe_load(VERDICTS_PATH.read_text())


@pytest.fixture(scope="module")
def campaign() -> dict[str, object]:
    return yaml.safe_load(CAMPAIGN_PATH.read_text())


@pytest.fixture(scope="module")
def prose() -> str:
    """The two deliverables that restate campaign numbers rather than rendering them."""
    docs = CAMPAIGN_PATH.resolve().parent.parent.parent
    return (docs / "report.tex").read_text() + (docs.parent / "README.md").read_text()


@pytest.fixture(scope="module")
def page(tmp_path_factory) -> str:
    """The page as actually built, from the committed tables the deliverables render from."""
    out = tmp_path_factory.mktemp("artifact") / "index.html"
    data = CAMPAIGN_PATH.resolve().parent.parent.parent.parent / "data"
    build_artifact.main(
        campaign=CAMPAIGN_PATH,
        verdicts=VERDICTS_PATH,
        pooled=data / "cofold_spread_pooled.csv",
        replicates=data / "cofold_replicates.csv",
        wildtype=data / "cofold_wildtype.csv",
        out=out,
    )
    return out.read_text()


class TestVerdictLedger:
    """One row per experiment, with a verdict the renderers know how to colour."""

    def test_every_experiment_has_the_fields_both_renderers_read(self, verdicts):
        for experiment in verdicts["experiments"]:
            missing = {"id", "what", "baseline", "gate", "result", "verdict"} - set(experiment)
            assert not missing, f"{experiment.get('id')} missing {missing}"

    def test_every_verdict_is_one_the_artifact_can_render(self, verdicts):
        """An unknown verdict renders as an unstyled pill rather than failing, so it has to be
        caught here or it ships looking like a bug."""
        found = {e["verdict"] for e in verdicts["experiments"]}
        assert found <= VERDICT_VALUES, f"unknown verdicts: {found - VERDICT_VALUES}"

    def test_the_not_established_list_is_not_empty(self, verdicts):
        """A campaign that claims nothing is unestablished has stopped looking."""
        assert len(verdicts["not_established"]) >= 3


class TestCampaignNumbers:
    """The measured values that are not a CSV row."""

    def test_each_stage_accounts_for_what_it_removed(self, campaign):
        """`surviving + removed` of one stage must be the previous stage's survivors, or the
        page shows a funnel whose bars do not add up. With `removed >= 0` this also makes the
        funnel monotonic, so nothing separately asserts that."""
        stages = campaign["funnel"]
        assert all(stage["removed"] >= 0 for stage in stages)
        for before, after in zip(stages, stages[1:]):
            if after["stage"] == "submitted":  # capped at 20, so it is not a filter
                continue
            assert after["surviving"] + after["removed"] == before["surviving"], after["stage"]

    def test_the_escape_cone_shrinks_with_reach(self, campaign):
        grid = campaign["escape_cone"]["by_reach_angstrom"]
        values = [grid[k] for k in sorted(grid)]
        assert values == sorted(values, reverse=True)

    def test_the_headline_agrees_with_the_funnel(self, campaign):
        """The strip at the top of the page and the funnel below it are the same numbers."""
        head, stages = campaign["headline"], {s["stage"]: s for s in campaign["funnel"]}
        assert head["screened"] == stages["engagement"]["surviving"]
        assert head["novelty_admissible"] == stages["fold novelty"]["surviving"]
        assert head["above_bar"] + head["tied_bar"] == stages["at or above the receptor's bar"]["surviving"]
        assert head["clear_zero"] == stages["clears zero on the drop"]["surviving"]
        assert head["submitted"] == stages["submitted"]["surviving"]
        assert head["above_bar"] + head["tied_bar"] + head["below_bar"] == head["screened"]

    def test_tier_one_is_recorded_as_inverted(self, campaign):
        """The campaign's central negative result. If this ever goes negative the strip, the
        report and the ledger all have to be rewritten, so it fails loudly rather than quietly
        rendering a claim nobody re-read."""
        assert campaign["headline"]["tier1_control_bias"] > 0


class TestWandbReport:
    """The narrative half of the report, which survives every refresh and so can go stale."""

    @pytest.mark.parametrize(
        "builder",
        [
            build_wandb_report.highlights,
            build_wandb_report.task,
            build_wandb_report.learned,
            build_wandb_report.nxt,
        ],
    )
    def test_every_narrative_section_renders_from_the_campaign(self, builder, campaign) -> None:
        """Each reads campaign.yaml by key, so a renamed field is a KeyError here rather than
        at report time, when the only signal is a half-written report."""
        assert builder(campaign).startswith("## ")

    def test_the_run_log_carries_every_experiment(self, verdicts) -> None:
        """The collapsed log is the only place a DISCARD row is guaranteed to appear; a filter
        that drops one leaves the report claiming a cleaner campaign than was run."""
        log = build_wandb_report.run_log_markdown(verdicts)
        for experiment in verdicts["experiments"]:
            assert f"`{experiment['id']}`" in log

    def test_the_novelty_count_is_summed_from_the_funnel(self, campaign) -> None:
        """Two sections quote this count; both call the helper, so neither can drift."""
        removed = sum(
            s["removed"] for s in campaign["funnel"] if s["stage"].endswith("novelty")
        )
        assert build_wandb_report.novelty_removed(campaign) == removed
        assert f"removed {removed} of" in build_wandb_report.highlights(campaign)
        assert f"| {removed} of" in build_wandb_report.learned(campaign)


class TestArtifactPage:
    """The built page, which is a graded deliverable and had shipped with the control in it."""

    def test_the_variant_order_still_matches_the_measurement_stack(self) -> None:
        """The page keeps its own copy so it does not import gemmi/Biopython to render text. A
        variant added to `common.VARIANTS` and not here does not error -- its rows silently
        vanish from both tier-1 figures -- so the drift has to fail somewhere."""
        assert build_artifact.VARIANT_ORDER == tuple(common.VARIANTS)

    def test_the_control_is_never_drawn_as_a_design(self, page) -> None:
        """It is the bar every chart is measured against, so drawing it as one of the designs
        it is the bar for both miscounts the cohort and shows the excluded run as a result."""
        assert build_artifact.CONTROL_DESIGN not in page

    def test_the_page_carries_every_funnel_stage(self, campaign, page) -> None:
        for stage in campaign["funnel"]:
            assert build_artifact.esc(stage["stage"]) in page

    def test_the_page_carries_the_retraction(self, campaign, page) -> None:
        """The unpaired pass's claim and its paired refutation both have to stay on the page:
        a deliverable that quietly drops a withdrawn result is the failure this campaign is
        about."""
        claim = campaign["tier1_paired"]["withdrawn_claim"]["bulge_core_fisher_p"]
        result = campaign["tier1_paired"]["paired_result"]["bulge_core_fisher_p"]
        assert f"{claim:.4f}" in page and f"{result:.4f}" in page

    def test_the_scoreboard_draws_one_bar_per_screened_design(self, campaign, page) -> None:
        """The caption states the denominator, so the bars have to match it."""
        chart = page.split("screened design against the receptor bar")[1].split("</svg>")[0]
        assert chart.count("<rect") == campaign["headline"]["screened"]


class TestHandWrittenDeliverables:
    """The README and the LaTeX report restate campaign numbers instead of rendering them.

    Generating them from the YAML would be the deeper fix, but the report is written once and
    submitted, so it would buy a build risk at the deadline against a drift risk. Detection is
    what is actually needed: every load-bearing number has to still appear, as written, in both.
    """

    @pytest.mark.parametrize(
        ("path", "fmt"),
        [
            ("headline.bar", "{:.3f}"),
            ("headline.best_drop", "{:.3f}"),
            ("headline.tier1_control_bias", "{:.4f}"),
            ("headline.tier1_control_bias_sem", "{:.4f}"),
            ("headline.screened", "{}"),
            ("headline.clear_zero", "{}"),
            ("headline.tied_bar", "{}"),
            ("escape_cone.free_rays_at_8A", "{}"),
            ("variance.designs.between_sd", "{}"),
            ("variance.designs.within_sem", "{}"),
            ("bar.invocations", "{}"),
            ("tier1_paired.control_bar", "{:.4f}"),
            ("tier1_paired.cohort_mean", "{:.4f}"),
            ("tier1_paired.predictions", "{}"),
            ("tier1_paired.paired_result.bulge_core_fisher_p", "{:.4f}"),
            ("tier1_paired.withdrawn_claim.bulge_core_fisher_p", "{:.4f}"),
        ],
    )
    def test_every_load_bearing_number_still_appears_as_written(self, path, fmt, campaign, prose):
        node = campaign
        for part in path.split("."):
            node = node[part]
        assert fmt.format(node) in prose, f"{path} = {node} is no longer in the report or README"
