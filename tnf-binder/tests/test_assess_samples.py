"""Reporting whether a structural verdict survives the sample it was computed on.

The clash counts are the real per-sample numbers from the two slot-filling rounds
(`docs/clashes-per-sample.txt`), so the suite fails if this stops reproducing the 10-of-16
flip rate that prompted it.
"""

from __future__ import annotations


import pytest

from tnf_binder.assess import (
    Assessment,
    SampledAssessment,
    geometry_across_samples,
    unstable_verdicts,
)
from tnf_binder.verifiers.interface import InterfaceMetrics

# design -> clashes on model_0, model_1, model_2, measured on all three structures
PER_SAMPLE = {
    "helix0_E67H_E70H_R72H": (1, 1, 4),
    "helix0_R21H_E70H_R72H": (7, 3, 7),
    "helix0_R21H_R35H_E67H": (1, 0, 1),
    "helix0_R21H_R35H_E70H": (3, 2, 0),
    "helix0_R21H_R35H_R72H": (0, 1, 4),
    "helix0_R35H_E67H_E70H": (0, 1, 9),
    "helix0_R35H_E67H_R72H": (0, 0, 4),
    "helix0_R35H_E70H_R72H": (1, 0, 6),
}
# interface area for the design the submission took on the strength of its model_0
AREAS = (658.2, 837.6, 1203.2)


def sampled(name: str, clashes, areas=None) -> SampledAssessment:
    areas = areas or [600.0] * len(clashes)
    return SampledAssessment(name, tuple(
        Assessment(name, metrics=InterfaceMetrics(
            area=a, epitope_focus=0.5, divergent_fraction=0.2, clashes=c))
        for c, a in zip(clashes, areas)))


class TestSampledAssessment:
    """One design across its samples, reported rather than re-judged."""

    def test_a_stable_verdict_has_no_caveat(self):
        assert sampled("x", (1, 1, 4)).caveat() is None
        assert sampled("x", (0, 0, 0)).caveat() is None

    def test_the_submitted_design_had_a_sample_dependent_verdict(self):
        """tnf2_13 was taken on model_0 reading 0 clashes; model_2 reads 9."""
        design = sampled("helix0_R35H_E67H_E70H", (0, 1, 9), AREAS)
        assert not design.clash_free_is_stable
        caveat = design.caveat()
        assert "sample-dependent" in caveat and "[0, 1, 9]" in caveat

    def test_the_caveat_carries_the_area_spread(self):
        """Area moved 658 to 1203 A^2 on the same design — a factor the verdict hid."""
        assert "545" in sampled("x", (0, 1, 9), AREAS).caveat()

    def test_clashes_are_named_for_the_design_so_a_report_can_identify_it(self):
        assert sampled("helix0_R72H", (0, 1)).clashes.name == "helix0_R72H"


class TestUnstableVerdicts:
    """The line to print before quoting a count from a structural gate."""

    def test_reproduces_the_campaign_flip_rate(self):
        report = unstable_verdicts([sampled(n, c) for n, c in PER_SAMPLE.items()])
        assert "6 sample-dependent (75%)" in report and "NOT trustworthy" in report

    def test_names_every_design_whose_verdict_would_change(self):
        report = unstable_verdicts([sampled(n, c) for n, c in PER_SAMPLE.items()])
        for flipped in ("helix0_R21H_R35H_E67H", "helix0_R35H_E67H_E70H"):
            assert flipped in report
        assert "helix0_R21H_E70H_R72H" not in report  # 7,3,7 — never clash-free, stable

    def test_a_wholly_stable_set_says_so_without_caveats(self):
        report = unstable_verdicts([sampled("a", (0, 0, 0)), sampled("b", (2, 2, 2))])
        assert "0 sample-dependent" in report and "NOT trustworthy" not in report


class TestGeometryAcrossSamples:
    """Merging score_designs rows per design, at the least favourable sample.

    Numbers are the real ones for `helix0_R21H_E38H_E67H` from the A38 round, whose clash-free
    verdict and epitope focus both depend on which structure is read.
    """

    ROWS = [
        {"design": "A__helix0_R21H_E38H_E67H_model_0", "clashes": "0",
         "interface_area": "629.7", "epitope_focus": "0.48", "flag": ""},
        {"design": "A__helix0_R21H_E38H_E67H_model_1", "clashes": "0",
         "interface_area": "1333.0", "epitope_focus": "0.24", "flag": "off_epitope"},
        {"design": "A__helix0_R21H_E38H_E67H_model_2", "clashes": "1",
         "interface_area": "1303.0", "epitope_focus": "0.22", "flag": "clash(1),off_epitope"},
    ]

    def test_merges_every_sample_into_one_row_per_design(self):
        merged = geometry_across_samples(self.ROWS)
        assert list(merged) == ["helix0_R21H_E38H_E67H"]
        assert merged["helix0_R21H_E38H_E67H"]["n_samples"] == "3"

    def test_takes_the_least_favourable_sample_of_each_metric(self):
        """Worst case, the same convention delta_q_floor uses: most clashes, least area,
        least focus. Reading model_0 alone would have called this one clash-free at 0.48."""
        g = geometry_across_samples(self.ROWS)["helix0_R21H_E38H_E67H"]
        assert g["clashes"] == "1"
        assert float(g["interface_area"]) == pytest.approx(629.7)
        assert float(g["epitope_focus"]) == pytest.approx(0.22)

    def test_a_flag_raised_on_any_sample_is_carried(self):
        g = geometry_across_samples(self.ROWS)["helix0_R21H_E38H_E67H"]
        assert "off_epitope" in g["flag"] and "clash(1)" in g["flag"]

    def test_records_the_spread_so_a_caller_can_say_the_verdict_was_close(self):
        g = geometry_across_samples(self.ROWS)["helix0_R21H_E38H_E67H"]
        assert g["sample_spread"] == "0-1" and g["focus_spread"] == "0.22-0.48"

    def test_a_design_with_one_sample_is_unchanged_and_marked_as_one(self):
        g = geometry_across_samples(self.ROWS[:1])["helix0_R21H_E38H_E67H"]
        assert g["n_samples"] == "1" and g["clashes"] == "0" and g["sample_spread"] == "0-0"
