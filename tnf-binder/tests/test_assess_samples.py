"""Reporting whether a structural verdict survives the sample it was computed on.

The clash counts are the real per-sample numbers from the two slot-filling rounds
(`docs/clashes-per-sample.txt`), so the suite fails if this stops reproducing the 10-of-16
flip rate that prompted it.
"""

from __future__ import annotations


from tnf_binder.assess import Assessment, SampledAssessment, unstable_verdicts
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
