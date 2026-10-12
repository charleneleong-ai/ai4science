"""Whether a verdict survives the sample it was computed on.

The numbers are the real per-sample clash counts from the campaign this was written after, so
a change in the logic contradicts something observed rather than failing a synthetic case.
"""

from __future__ import annotations

import pytest

from touchstone.sampled import Sampled, stability

# clashes per diffusion sample, scored on all three structures of each design
CLASHES = {
    "E67H_E70H_R72H": (1, 1, 4),      # never clash-free: stable drop
    "R21H_E70H_R72H": (7, 3, 7),      # stable drop
    "R21H_R35H_E67H": (1, 0, 1),      # FLIPS
    "R21H_R35H_E70H": (3, 2, 0),      # FLIPS
    "R21H_R35H_R72H": (0, 1, 4),      # FLIPS
    "R35H_E67H_E70H": (0, 1, 9),      # FLIPS — submitted on the strength of model_0
    "R35H_E67H_R72H": (0, 0, 4),      # FLIPS
    "R35H_E70H_R72H": (1, 0, 6),      # FLIPS
}
MEASURED = [Sampled.of(name, values) for name, values in CLASHES.items()]


def clash_free(value: float) -> bool:
    return value == 0


class TestSampled:
    """One quantity across the samples of one design."""

    def test_worst_is_the_least_favourable_sample(self):
        """A conservative gate reads this, not the mean of a count nothing measures."""
        design = Sampled.of("R35H_E67H_E70H", (0, 1, 9))
        assert design.worst() == 9
        assert design.worst(higher_is_worse=False) == 0

    def test_spread_exposes_what_one_sample_hides(self):
        assert Sampled.of("x", (0, 1, 9)).spread == 9
        assert Sampled.of("x", (1, 1, 1)).spread == 0

    @pytest.mark.parametrize("values,verdict", [
        ((0, 0, 0), "keep"),
        ((1, 1, 4), "drop"),
        ((0, 1, 9), "unstable"),
        ((), "unmeasured"),
    ])
    def test_verdict_names_the_coin_flip(self, values, verdict):
        assert Sampled.of("x", values).verdict(clash_free) == verdict

    def test_agreement_is_the_fraction_that_say_keep(self):
        assert Sampled.of("x", (0, 0, 9)).agreement(clash_free) == pytest.approx(2 / 3)

    def test_an_unmeasured_design_is_not_a_passing_one(self):
        empty = Sampled.of("x", ())
        assert not empty and empty.agreement(clash_free) == 0.0 and empty.mean == 0.0


class TestStability:
    """Whether a gate's output over a set of designs says anything about the designs."""

    def test_counts_the_real_campaign_flips(self):
        """6 of these 8 flip — which is why "six of eight failed on structure" was a statement
        about model_0 rather than about the designs."""
        report = stability(MEASURED, clash_free)
        assert len(report.unstable) == 6
        assert set(report.drop) == {"E67H_E70H_R72H", "R21H_E70H_R72H"}
        assert report.keep == ()

    def test_a_gate_deciding_most_designs_by_coin_flip_is_not_trustworthy(self):
        report = stability(MEASURED, clash_free)
        assert report.unstable_fraction == pytest.approx(0.75)
        assert not report.trustworthy and "NOT trustworthy" in str(report)

    def test_a_stable_gate_is_trustworthy(self):
        """Same designs, a threshold loose enough that no sample disagrees."""
        report = stability(MEASURED, lambda v: v <= 9)
        assert report.unstable == () and report.trustworthy

    def test_every_design_is_accounted_for(self):
        report = stability(MEASURED, clash_free)
        assert report.total == len(MEASURED)
