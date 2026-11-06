"""Selection logic and co-fold readout — the parts that decide which designs get submitted."""

from __future__ import annotations

import json
import math

import pytest
import typer
from cofold_followup import (
    SIGMA,
    Contrast,
    Paired,
    Spread,
    WildtypeRow,
    comparison,
    COHORT_BASIS,
    Cohort,
    Pooled,
    difference_sem,
    invocation_predictions,
    paired_drops,
    pooled_between_sd,
    run_effect_sd,
    spanning,
    spread_of,
)
from cofold_followup import RunRow
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
    """A spread table as `report-spread` writes it, keyed by design, from (drop_mean, drop_sem).

    `drop_sd` is deliberately not equal to `drop_sem` -- they were once the same value here, and
    while they were, nothing could tell which column the consumer actually divided by.
    """
    return {
        design: {
            "design": design,
            "drop_mean": str(mean),
            "drop_sd": str(round(sem * math.sqrt(5), 6)),
            "drop_sem": str(sem),
            "sem_basis": "own_samples",
        }
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


class TestReplicateStatistics:
    """The between-run component, measured by separate invocations rather than more samples."""

    def test_the_between_run_spread_is_of_the_run_means_not_of_the_samples(self):
        """The eight invocations of the control, whose spread is what the pass exists to get."""
        bar = Paired("control_tcr", (0.100, 0.095, 0.104, 0.099, 0.108, 0.099, 0.102, 0.098))
        assert bar.mean == pytest.approx(0.10063, abs=1e-5)
        assert bar.between_sd == pytest.approx(0.00400, abs=1e-5)
        assert bar.sem == pytest.approx(bar.between_sd / math.sqrt(8))

    @pytest.mark.parametrize(
        "values,verdict",
        [
            ((0.100, 0.102, 0.098, 0.101), "above"),
            ((0.001, -0.002, 0.003, -0.001), "tied"),
            ((-0.05, -0.04, -0.06, -0.05), "below"),
        ],
    )
    def test_a_repeated_quantity_is_tested_against_zero_by_the_shared_rule(self, values, verdict):
        assert Paired("d", values).verdict == verdict

    @pytest.mark.parametrize(
        "between,within,expected",
        [
            # The control, measured: a spread smaller than its own sampling error leaves no
            # run effect to charge a design for, and must not come back as a hypot of the two.
            (0.004, 0.007, 0.0),
            (0.007, 0.007, 0.0),
            (0.013, 0.005, 0.012),
            (0.010, 0.000, 0.010),
        ],
    )
    def test_the_two_components_subtract_because_the_spread_already_holds_the_sampling(
        self, between, within, expected
    ):
        """Var(run means) = Var(run effect) + within**2, so adding them would double-count."""
        assert run_effect_sd(between, within) == pytest.approx(expected, abs=5e-4)

    def test_variances_pool_across_designs_not_standard_deviations(self):
        """Averaging SDs would understate a cohort whose designs vary in reproducibility."""
        sds = (0.02, 0.06)
        assert pooled_between_sd(sds) == pytest.approx(0.044721, abs=1e-6)
        assert pooled_between_sd(sds) > sum(sds) / 2


class TestInvocationPredictions:
    """Both shapes a replicated pass arrives in, so no caller has to know which it has."""

    def test_a_directory_of_invocations_resolves_each_ones_predictions(self, tmp_path):
        for name in ("rep2", "rep1"):
            (tmp_path / name / "boltz_results_x" / "predictions").mkdir(parents=True)
        found = invocation_predictions(tmp_path)
        assert [p.parent.parent.name for p in found] == ["rep1", "rep2"]
        assert all(p.name == "predictions" for p in found)

    def test_one_pass_out_dir_is_a_single_invocation(self, tmp_path):
        (tmp_path / "boltz_results_x" / "predictions").mkdir(parents=True)
        assert invocation_predictions(tmp_path) == [tmp_path / "boltz_results_x" / "predictions"]

    def test_a_sibling_holding_no_predictions_is_not_an_invocation(self, tmp_path):
        """Counting a stray directory would inflate the denominator of every between-run claim."""
        (tmp_path / "rep1" / "boltz_results_x" / "predictions").mkdir(parents=True)
        (tmp_path / "logs").mkdir()
        assert len(invocation_predictions(tmp_path)) == 1

    def test_refuses_a_directory_holding_neither(self, tmp_path):
        with pytest.raises(typer.BadParameter, match="no predictions"):
            invocation_predictions(tmp_path)


class TestCohortDenominator:
    """Why a two-run design borrows the field's spread instead of trusting its own two points."""

    @pytest.mark.parametrize(
        "drops,expected",
        [
            # Two runs that happen to agree report zero spread; taken at face value that makes
            # any drop infinitely significant, which is how +0.011 came to read as "above".
            ((0.011, 0.011), "tied"),
            ((0.112, 0.112), "above"),
            ((-0.093, -0.093), "below"),
        ],
    )
    def test_a_designs_own_spread_is_not_its_error_at_two_runs(self, drops, expected):
        cohort = Cohort(between_sd=0.018, within_sem=0.015)
        assert Paired("d", drops).between_sd == pytest.approx(0.0)
        assert cohort.verdict(Paired("d", drops)) == expected

    def test_the_pooled_sem_tightens_with_the_number_of_invocations(self):
        cohort = Cohort(between_sd=0.018, within_sem=0.015)
        assert cohort.sem(2) == pytest.approx(0.018 / math.sqrt(2))
        assert cohort.sem(8) < cohort.sem(2)

    def test_the_run_effect_is_what_the_spread_leaves_over_the_sampling(self):
        assert Cohort(between_sd=0.018, within_sem=0.015).run_effect == pytest.approx(
            math.sqrt(0.018**2 - 0.015**2), abs=1e-6
        )

    def test_the_control_is_excluded_from_the_spread_it_is_the_bar_for(self):
        """Pooling the bar's own irreproducibility in would let it widen its own margin."""
        def runs(design, drops):
            return [
                RunRow(replicate=f"r{i}", spread=Spread(design, (d + 0.7,) * 5, (0.7,) * 5), wildtype_mean=None)
                for i, d in enumerate(drops)
            ]

        by_design = {CONTROL: runs(CONTROL, (0.10, 0.15)), "d": runs("d", (0.05, 0.07))}
        paired = paired_drops(by_design)
        cohort = Cohort.of(
            {d: [r.spread.drop_sem for r in rs] for d, rs in by_design.items()}, paired
        )
        # The control's two runs differ by 0.05, the design's by 0.02; only the design counts.
        assert cohort.between_sd == pytest.approx(paired["d"].between_sd, abs=1e-9)
        assert cohort.between_sd < paired[CONTROL].between_sd


class TestPooledTable:
    """The pooled table the submission ranks on, and the scale it declares in its own rows."""

    @staticmethod
    def pooled(drops: tuple[float, ...], between: float = 0.018) -> Pooled:
        return Pooled("d", Paired("d", drops), Cohort(between_sd=between, within_sem=0.015))

    def test_the_spread_column_is_the_cohorts_not_the_designs(self):
        """Two runs that agree report 0.000; ranking on that is the failure Cohort prevents."""
        row = self.pooled((0.112, 0.112)).row
        assert row["drop_sd"] == pytest.approx(0.018)
        assert row["drop_own_sd"] == pytest.approx(0.0)
        assert row["drop_sd"] != row["drop_own_sd"]

    def test_the_invariant_the_consumer_relies_on_still_holds(self):
        """`submission` reads drop_sd and drop_sem separately, so sem must stay sd/sqrt(n)."""
        row = self.pooled((0.10, 0.12)).row
        assert row["drop_sem"] == pytest.approx(row["drop_sd"] / math.sqrt(2), abs=1e-4)

    def test_every_row_names_the_scale_it_was_measured_on(self):
        """Otherwise a reader of the published CSV has to know which command wrote it."""
        assert self.pooled((0.10, 0.12)).row["sem_basis"] == COHORT_BASIS
        assert Spread(design="d", complexed=(0.9, 0.8), apo=(0.7, 0.6)).row["sem_basis"] != COHORT_BASIS

    def test_both_producers_write_the_keys_the_submission_reads(self):
        """The one contract candidates() depends on, asserted for both tables that can feed it."""
        needed = {"design", "drop_mean", "drop_sd", "drop_sem"}
        assert needed <= set(self.pooled((0.10, 0.12)).row)
        assert needed <= set(Spread(design="d", complexed=(0.9, 0.8), apo=(0.7, 0.6)).row)

    def test_the_two_routes_to_significance_agree(self):
        """`Pooled` decides via Cohort and `Candidate` recomputes from the columns; the module
        docstring promises the printed table and the uploaded file cannot disagree, so that must
        hold by test rather than by arithmetic coincidence."""
        row = self.pooled((0.112, 0.112)).row
        assert row["drop_significant"] is one_candidate(
            float(row["drop_mean"]), float(row["drop_sem"])
        ).drop_significant


class TestTier1Difference:
    """The wild-type difference: paired inside an invocation, scaled by the cohort."""

    @staticmethod
    def contrast(design: str, replicate: str, difference: float, spread: float = 0.002) -> Contrast:
        low = tuple(0.930 + spread * i for i in range(5))
        return Contrast(
            design=design,
            replicate=replicate,
            high=tuple(v + difference for v in low),
            low=low,
        )

    def test_pairing_inside_the_invocation_cancels_the_run_offset(self):
        """The defect this pass exists to fix. Two invocations that differ by a large offset in
        both conditions must give the same difference as two that do not -- an unpaired
        subtraction would report the offset as signal."""
        level = self.contrast("d", "r1", difference=-0.010)
        shifted = Contrast(
            design="d",
            replicate="r2",
            high=tuple(v + 0.050 for v in level.high),
            low=tuple(v + 0.050 for v in level.low),
        )
        assert shifted.mean == pytest.approx(level.mean, abs=1e-9)

    def test_the_error_is_the_two_conditions_in_quadrature(self):
        c = self.contrast("d", "r1", difference=-0.010)
        assert c.sem == pytest.approx(difference_sem(c.high, c.low), abs=1e-12)

    def test_the_drop_and_the_difference_share_one_quadrature_rule(self):
        """`Spread.drop_sem` and `Contrast.sem` are the same statistic; if they ever disagree the
        campaign has two denominators for one shape of quantity."""
        c = self.contrast("d", "r1", difference=-0.010)
        assert Spread(design="d", complexed=c.high, apo=c.low).drop_sem == pytest.approx(c.sem)

    @pytest.mark.parametrize(
        "design_diffs,expected",
        [
            ((-0.030, -0.031), "below"),  # clearly the other side of the control's bar
            ((+0.003, +0.002), "tied"),  # sitting on the bar
            ((+0.040, +0.041), "above"),  # scores wild type higher than even the control does
        ],
    )
    def test_a_design_is_judged_against_the_control_not_against_zero(self, design_diffs, expected):
        """The control's own difference is positive, so zero is the wrong reference: a small
        negative difference is not selectivity, it is less bias than the control shows."""
        bar_diffs = (+0.0027, +0.0027)
        by_design = {
            CONTROL: [self.contrast(CONTROL, f"r{i}", d) for i, d in enumerate(bar_diffs)],
            "bulge_core_09": [
                self.contrast("bulge_core_09", f"r{i}", d) for i, d in enumerate(design_diffs)
            ],
        }
        paired = {k: Paired(k, tuple(c.mean for c in v)) for k, v in by_design.items()}
        cohort = Cohort.of({k: [c.sem for c in v] for k, v in by_design.items()}, paired)
        row = WildtypeRow(
            design="bulge_core_09",
            differences=paired["bulge_core_09"],
            cohort=cohort,
            bar=paired[CONTROL],
        )
        assert row.versus_bar == expected

    def test_the_published_row_carries_the_bar_it_was_judged_against(self):
        """`versus_bar` is not reproducible from the file unless the bar's mean, its error and
        the number of invocations behind it all travel with the row."""
        by_design = {
            CONTROL: [self.contrast(CONTROL, f"r{i}", d) for i, d in enumerate((0.003, 0.002))],
            "hub_01": [self.contrast("hub_01", f"r{i}", d) for i, d in enumerate((0.030, 0.031))],
        }
        paired = {k: Paired(k, tuple(c.mean for c in v)) for k, v in by_design.items()}
        cohort = Cohort.of({k: [c.sem for c in v] for k, v in by_design.items()}, paired)
        row = WildtypeRow(
            design="hub_01", differences=paired["hub_01"], cohort=cohort, bar=paired[CONTROL]
        ).row
        assert {"control_bar", "control_bar_sem", "control_bar_runs"} <= set(row)
        assert row["control_bar_runs"] == 2
        assert row["variant"] == "hub"
        assert row["sem_basis"] == COHORT_BASIS

    def test_the_cohort_denominator_is_condition_agnostic(self):
        """`Cohort.of` takes the within-run errors as numbers, so the drop and the tier-1
        difference borrow one pooling rule instead of each growing its own."""
        paired = {"a": Paired("a", (0.10, 0.12)), "b": Paired("b", (0.20, 0.26))}
        cohort = Cohort.of({"a": [0.01, 0.01], "b": [0.02, 0.02]}, paired)
        assert cohort.between_sd == pytest.approx(
            pooled_between_sd([paired["a"].between_sd, paired["b"].between_sd])
        )
        assert cohort.within_sem == pytest.approx(0.015)
