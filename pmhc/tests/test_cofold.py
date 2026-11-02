"""Selection logic and co-fold readout — the parts that decide which designs get submitted."""

from __future__ import annotations

import json

import pytest
from cofold_scores import Dependence, collect, peptide_pair
from cofold_specs import (
    CONTROL_PEPTIDE_SHARE,
    FOLD_SEVERITY,
    chain_ids,
    conditions_from,
    judge,
    parse_run_name,
    run_name,
    spec_text,
)
from fold_novelty import FOLD_VERDICTS


def novelty_row(novel: bool = True, copied: float = 0.05) -> dict[str, str]:
    return {"novel": str(novel), "copied_fraction": str(copied)}


def engagement_row(share: float = 0.3, touched: int = 3, nominated: int = 3) -> dict[str, str]:
    return {
        "variant": "bulge_core",
        "peptide_share": str(share),
        "hotspots_touched": str(touched),
        "hotspots_nominated": str(nominated),
        "mutation_contacts": "6",
    }


class TestJudge:
    @pytest.mark.parametrize(
        "nov,fold,eng,expected",
        [
            (novelty_row(), "shared_fold", engagement_row(), "selected"),
            (novelty_row(False, 0.9), "shared_fold", engagement_row(), "rejected_sequence"),
            (novelty_row(), "copy_of_known", engagement_row(), "rejected_copy_of_known"),
            (novelty_row(), "ig_sandwich", engagement_row(), "rejected_ig_sandwich"),
            (novelty_row(), "shared_fold", engagement_row(share=0.05), "rejected_groove_gripper"),
            (novelty_row(), "shared_fold", engagement_row(touched=2), "rejected_missed_hotspot"),
        ],
    )
    def test_each_gate_has_its_own_verdict(self, nov, fold, eng, expected):
        """The reason is the point: a funnel that only counts survivors cannot be audited."""
        assert judge("d", nov, fold, eng).verdict == expected

    def test_the_most_serious_reason_wins(self):
        """A plagiarised groove-gripper is reported as plagiarised, not as a groove-gripper."""
        fate = judge("d", novelty_row(False, 0.9), "copy_of_known", engagement_row(share=0.0))
        assert fate.verdict == "rejected_sequence"

    def test_the_share_bar_is_the_measured_receptor(self):
        """Not a pinned literal: the bar is 7OW6's own split, 55 peptide contacts against 220
        on the HLA, so it cannot drift away from the measurement it came from."""
        assert CONTROL_PEPTIDE_SHARE == pytest.approx(55 / (55 + 220))
        assert judge("d", novelty_row(), "shared_fold", engagement_row(share=0.200)).verdict == (
            "selected"
        )
        assert judge("d", novelty_row(), "shared_fold", engagement_row(share=0.199)).verdict == (
            "rejected_groove_gripper"
        )


class TestFoldSeverity:
    def test_covers_every_verdict_the_producer_can_return(self):
        """worst_fold indexes this by verdict string, so a verdict missing here is a KeyError
        mid-selection rather than a wrong answer."""
        assert set(FOLD_SEVERITY) == set(FOLD_VERDICTS)

    def test_orders_a_copy_above_a_shared_fold(self):
        assert FOLD_SEVERITY["copy_of_known"] > FOLD_SEVERITY["ig_sandwich"]
        assert FOLD_SEVERITY["ig_sandwich"] > FOLD_SEVERITY["shared_fold"]
        assert FOLD_SEVERITY["shared_fold"] > FOLD_SEVERITY["distinct_fold"]


class TestRunNames:
    @pytest.mark.parametrize(
        "design,condition",
        [("bulge_core_13", "complex_g12d"), ("control_tcr", "apo_hla")],
    )
    def test_round_trips(self, design, condition):
        assert parse_run_name(run_name(design, condition)) == (design, condition)

    def test_a_name_without_a_separator_is_not_a_run(self):
        assert parse_run_name("gen_out") is None


class TestChainIds:
    def test_a_design_is_one_chain_before_the_condition(self):
        assert chain_ids("bulge_core_13", ("A", "B", "C")) == ("X", "A", "B", "C")

    def test_the_control_is_two(self):
        assert chain_ids("control_tcr", ("A", "B", "C")) == ("D", "E", "A", "B", "C")

    def test_apo_has_no_peptide_chain(self):
        assert "C" not in chain_ids("bulge_core_13", ("A", "B"))

    def test_reads_the_condition_chains_in_spec_order(self, tmp_path):
        targets = tmp_path / "cofold_targets.json"
        targets.write_text(
            json.dumps(
                [
                    {"name": "complex_g12d", "chains": {"A": "GG", "B": "KK", "C": "VV"}},
                    {"name": "apo_hla", "chains": {"A": "GG", "B": "KK"}},
                ]
            )
        )
        chains = conditions_from(targets)
        assert chain_ids("d", tuple(chains["complex_g12d"])) == ("X", "A", "B", "C")
        assert chain_ids("d", tuple(chains["apo_hla"])) == ("X", "A", "B")


class TestSpecText:
    def test_every_chain_gets_an_empty_msa(self):
        """--use_msa_server is a flag, not a value; a spec without this silently seeks an MSA."""
        assert spec_text([("X", "AAA"), ("A", "GGG")]).count("msa: empty") == 2


class TestPeptidePair:
    def test_averages_the_binder_rows_against_the_peptide_column(self):
        confidence = {"pair_chains_iptm": {"0": {"3": 0.90}, "1": {"3": 0.80}, "2": {"3": 0.10}}}
        assert peptide_pair(confidence, ("D", "E", "A", "C"), ("D", "E")) == pytest.approx(0.85)

    def test_averages_both_directions_of_an_asymmetric_matrix(self):
        """Boltz's pair matrix is not symmetric and neither direction is privileged."""
        confidence = {"pair_chains_iptm": {"0": {"3": 0.90}, "3": {"0": 0.70}}}
        assert peptide_pair(confidence, ("X", "A", "B", "C"), ("X",)) == pytest.approx(0.80)

    @pytest.mark.parametrize(
        "confidence,order",
        [
            # The apo condition has no peptide, so there is nothing to pair against.
            ({"pair_chains_iptm": {"0": {"1": 0.9}}}, ("X", "A", "B")),
            ({}, ("X", "A", "C")),
        ],
    )
    def test_none_when_there_is_nothing_to_report(self, confidence, order):
        assert peptide_pair(confidence, order, ("X",)) is None


class TestDependence:
    @pytest.mark.parametrize(
        "complexed,apo,expected",
        [(0.80, 0.30, 0.50), (0.80, 0.81, -0.01)],
    )
    def test_drop_is_the_confidence_lost_with_the_peptide_removed(self, complexed, apo, expected):
        row = Dependence("d", complexed, apo, None, None, 0.9)
        assert row.drop == pytest.approx(expected)
        assert row.row["peptide_dependence_drop"] == pytest.approx(expected)

    def test_absent_optional_conditions_stay_null_rather_than_zero(self):
        """A missing wild-type run must not read as a wild-type ipTM of 0."""
        row = Dependence("d", 0.8, 0.3, None, None, 0.9).row
        assert row["iptm_complex_wt"] is None
        assert row["design_peptide_pair_iptm"] is None


class TestCollect:
    def test_counts_every_directory_including_the_ones_it_cannot_score(self, tmp_path):
        """Boltz leaves a directory behind for a failed example. It must not become a score,
        and it must not vanish from the denominator either."""
        preds = tmp_path / "preds"
        preds.mkdir()
        for name in ("d__complex_g12d", "d__apo_hla", "e__complex_g12d", "stray"):
            (preds / name).mkdir()
        for name in ("d__complex_g12d", "e__complex_g12d"):
            (preds / name / "confidence_x_model_0.json").write_text(json.dumps({"iptm": 0.7}))

        screen = collect(preds)
        assert set(screen.runs["d"]) == {"complex_g12d"}
        assert screen.directories == 4
        assert screen.unscored == 1
        assert screen.unnamed == 1
