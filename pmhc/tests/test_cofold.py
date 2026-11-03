"""Selection logic and co-fold readout — the parts that decide which designs get submitted."""

from __future__ import annotations

import json
import math

import pytest
import typer
from cofold_followup import SIGMA, Spread, comparison, spanning, spread_of
from cofold_scores import Dependence, collect, peptide_pair
from cofold_specs import (
    CONTROL,
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
import submission
from submission import Candidate, candidates, csv_problems, row_problem


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


def submission_row(**changes: str) -> dict[str, str]:
    return {"name": "pmhc_01", "sequence": "ACDEFGHIKL", "molecule_class": "single_chain"} | changes


def spread_rows(drops: dict[str, tuple[float, float]]) -> dict[str, dict[str, str]]:
    """A spread table as `report-spread` writes it, keyed by design, from (drop_mean, drop_sem)."""
    return {
        design: {"design": design, "drop_mean": str(mean), "drop_sd": str(sem), "drop_sem": str(sem)}
        for design, (mean, sem) in drops.items()
    }


def one_candidate(mean: float, sem: float, control: tuple[float, float] = (0.10, 0.01)) -> Candidate:
    return candidates(
        spread_rows({"a": (mean, sem), CONTROL: control}), {"a": "A" * 60}, {}, count=20
    )[0]


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


class TestSpanningSubset:
    RANKED = [(f"d{i}", 1.0 - 0.01 * i) for i in range(78)]

    @pytest.mark.parametrize("count,head", [(12, 4), (9, 0), (5, 4), (2, 1), (3, 5), (1, 0)])
    def test_covers_both_ends_of_the_ranking(self, count, head):
        """A noise estimate taken only from the top would not bound the ties further down.

        A head at least as large as the count is the case that used to break this silently,
        returning the top `count` and spanning nothing.
        """
        picked = spanning(self.RANKED, count, head)
        assert len(picked) == count
        assert picked[0] == "d0"
        assert picked[-1] == ("d77" if count > 1 else "d0")

    def test_keeps_the_contenders_whatever_the_spread(self):
        """The decision the estimate has to settle is whether the leaders are really ahead."""
        assert spanning(self.RANKED, 12, 4)[:4] == ["d0", "d1", "d2", "d3"]

    def test_asking_for_more_than_exists_returns_everything_in_order(self):
        assert spanning(self.RANKED, 200, 4) == [f"d{i}" for i in range(78)]


class TestSpread:
    def test_the_drop_spread_combines_both_conditions(self):
        """Each condition is sampled independently, so the difference is noisier than either.

        The number is pinned rather than recomputed from the implementation's own formula: an
        assertion that restates `hypot(sd, sd)` would pass just as happily if the combination
        rule were a plain sum.
        """
        spread = spread_of(
            "d", {"complex_g12d": (0.95, 0.93, 0.97), "apo_hla": (0.85, 0.82, 0.86)}
        )
        assert spread.drop_mean == pytest.approx(0.1067, abs=1e-4)
        assert spread.drop_sd == pytest.approx(0.0289, abs=1e-4)
        assert spread.drop_sd > max(spread.complex_sd, spread.apo_sd)

    def test_the_range_pairs_the_best_against_the_worst(self):
        """No complex draw belongs to any particular apo draw, so the extremes cross them."""
        spread = spread_of(
            "d", {"complex_g12d": (0.95, 0.93, 0.97), "apo_hla": (0.85, 0.82, 0.86)}
        )
        assert (spread.drop_min, spread.drop_max) == pytest.approx((0.93 - 0.86, 0.97 - 0.82))

    def test_reports_both_denominators(self):
        """One n would hide that the two conditions can come back with different sample counts."""
        spread = spread_of("d", {"complex_g12d": (0.9, 0.91, 0.92), "apo_hla": (0.8, 0.81)})
        assert spread.samples == "3/2"

    @pytest.mark.parametrize(
        "conditions",
        [
            {"complex_g12d": (0.9,), "apo_hla": (0.8,)},
            {"complex_g12d": (0.9, 0.91)},
            {},
        ],
    )
    def test_one_sample_has_no_spread_to_report(self, conditions):
        """Returning a 0.0 spread here would read as "measured, and tiny"."""
        assert spread_of("d", conditions) is None


class TestSubmissionFormat:
    def test_accepts_a_legal_row(self):
        """`row_problem` returns None for a clean row, which no reject case can assert."""
        assert row_problem(submission_row()) is None

    @pytest.mark.parametrize(
        "change,expected",
        [
            ({"molecule_class": "protein"}, "bad class"),
            ({"sequence": "ACDEFGHIKX"}, "non-standard residues"),
            ({"sequence": "ACDEFGHIK"}, "length 9"),
            ({"sequence": "A" * 251}, "length 251"),
        ],
    )
    def test_rejects_what_the_submit_page_rejects(self, change, expected):
        """`protein` is not an accepted class and 10-250 aa is the stated range."""
        assert expected in row_problem(submission_row(**change))

    @pytest.mark.parametrize(
        "rows,expected",
        [
            ([], "empty"),
            ([submission_row(), submission_row(name="pmhc_02")], "duplicate sequences"),
            ([submission_row(), submission_row(sequence="C" * 60)], "duplicate names"),
            ([submission_row(name=f"pmhc_{i:02d}", sequence="A" * (10 + i)) for i in range(21)],
             "21 rows > 20 limit"),
        ],
    )
    def test_file_level_rules_a_single_row_cannot_see(self, rows, expected):
        assert any(expected in p for p in csv_problems(rows))

    def test_a_wrong_header_is_reported_rather_than_raised(self):
        """`check` reads files this code did not write, and `row_problem` indexes the columns
        it checks -- so the column problem has to come back alone."""
        problems = csv_problems([{"id": "x", "seq": "ACDEFGHIKL", "kind": "single_chain"}])
        assert len(problems) == 1 and "columns" in problems[0]


class TestCandidateSelection:
    def test_ranks_by_mean_drop_and_names_in_rank_order(self):
        chosen = candidates(
            spread_rows({"a": (0.05, 0.01), "b": (0.12, 0.01), CONTROL: (0.10, 0.01)}),
            {"a": "A" * 60, "b": "C" * 60},
            {},
            count=20,
        )
        assert [(c.name, c.design) for c in chosen] == [("pmhc_01", "b"), ("pmhc_02", "a")]

    def test_submits_no_more_than_the_cap(self):
        """The cap defines the submission, so the truncation it implies has to be exercised."""
        drops = {f"d{i:02d}": (0.10 - i / 1000, 0.001) for i in range(25)}
        chosen = candidates(
            spread_rows({**drops, CONTROL: (0.10, 0.01)}),
            {d: "A" * (20 + i) for i, d in enumerate(drops)},
            {},
            count=3,
        )
        assert [c.design for c in chosen] == ["d00", "d01", "d02"]

    @pytest.mark.parametrize(
        "spread,sequences,message",
        [
            ({"a": (0.1, 0.01)}, {"a": "A" * 60}, CONTROL),
            ({CONTROL: (0.1, 0.01)}, {}, "control and no designs"),
            ({"a": (0.1, 0.01), CONTROL: (0.1, 0.01)}, {}, "no sequence"),
        ],
    )
    def test_refuses_to_rank_what_it_cannot_rank(self, spread, sequences, message):
        """Without a control the drops have no scale, which is the basis of the whole claim."""
        with pytest.raises(typer.BadParameter, match=message):
            candidates(spread_rows(spread), sequences, {}, count=20)

    @pytest.mark.parametrize(
        "mean,sd,significant",
        [(0.12, 0.01, True), (0.02, 0.08, False), (-0.05, 0.01, False)],
    )
    def test_a_drop_must_clear_the_standard_error_of_its_own_mean(self, mean, sd, significant):
        """The pilot's third and fourth place reversed sign at five samples; a drop smaller
        than its own noise is not a drop."""
        assert one_candidate(mean, sd).drop_significant is significant

    @pytest.mark.parametrize(
        "mean,sd,verdict",
        [(0.30, 0.01, "above"), (0.105, 0.01, "tied"), (0.00, 0.01, "below")],
    )
    def test_compared_against_the_control_as_a_band_not_a_line(self, mean, sd, verdict):
        assert one_candidate(mean, sd).versus_control == verdict

    def test_the_sigma_rule_is_the_one_the_spread_report_states(self):
        """Two modules computing it independently is how the printed table and the uploaded
        file come to disagree about which designs are tied."""
        assert submission.comparison is comparison and submission.SIGMA == SIGMA


class TestSpreadStatistics:
    def test_the_tested_quantity_is_a_mean_so_its_error_is_the_standard_error(self):
        """drop_sd is the spread of one difference; drop_mean is a mean of n, and testing a
        mean against the per-draw spread is a sqrt(n)-times-too-strict gate."""
        spread = Spread(design="d", complexed=(0.9, 0.8, 0.85, 0.82, 0.88), apo=(0.7,) * 5)
        assert spread.drop_sem == pytest.approx(spread.drop_sd / math.sqrt(5))

    def test_the_consumer_reads_the_keys_the_producer_writes(self):
        """`candidates()` reads these back off the CSV, so a rename here must break a test
        rather than only break production."""
        written = Spread(design="d", complexed=(0.9, 0.8), apo=(0.7, 0.6)).row
        assert {"design", "drop_mean", "drop_sd", "drop_sem"} <= set(written)
