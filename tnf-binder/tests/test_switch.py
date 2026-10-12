"""Building and measuring a pH switch, with both silent failures made separable.

A mutation that does not take and a residue PROPKA declines to score look identical downstream
if they are collapsed, and a whole scan once came back with every pKa missing because
histidines had been renamed without building the imidazole. PyRosetta and PROPKA are injected
here so the sequence is tested without either installed.

pKa values are the real ones measured on this campaign's lead chassis.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path

import pytest

from tnf_binder.switch import (
    MIDPOINT,
    HisMeasurement,
    HisVariant,
    SampledHis,
    SwitchProfile,
    delta_q,
    his_combinations,
    his_tag,
    measure_his,
    position_list,
    protonated_fraction,
    switch_profiles,
)

# measured on screen_pocket_helix_0: two ideal, one usable, one dead
MEASURED = {("A", 21): 7.05, ("A", 67): 7.11, ("A", 35): 6.07, ("A", 38): 5.35}

# Rows as pka_robustness.py wrote them for the three submitted designs that span the delta-Q
# range, with the delta_q_floor the submission reported for each. These are the real campaign
# numbers, so they pin the metric to what was actually submitted rather than to a re-derivation.
PUBLISHED = {
    "helix0_R21H_E67H_E70H": ([(21, 6.61, 7.55), (67, 6.74, 7.32), (70, 6.06, 6.14)], 1.38),
    "helix0_R72H": ([(72, 6.23, 6.50)], 0.57),
    "helix0_R21H": ([(21, 3.99, 6.70)], 0.01),
}


def robustness_rows(design: str) -> list[dict[str, str]]:
    return [{"design": design, "position": str(p), "pka_min": str(lo), "pka_max": str(hi)}
            for p, lo, hi in PUBLISHED[design][0]]


class FakeMutator:
    """Stands in for HisMutator: records what it was asked to do, never loads PyRosetta."""

    def __init__(self, *, succeeds: bool = True) -> None:
        self.succeeds = succeeds
        self.written: list[Path] = []

    def load(self, path: Path) -> str:
        return f"pose:{path.stem}"

    def mutate(self, pose: str, positions: list[tuple[str, int]]) -> bool:
        return self.succeeds

    def write(self, pose: str, dest: Path) -> None:
        dest.write_text("MUTANT\n")
        self.written.append(dest)


def fake_propka(_mutant: Path, _workdir: Path, *, table: dict) -> dict:
    return table


class TestNaming:
    """The tag three scripts had each spelled out."""

    @pytest.mark.parametrize("positions,tag,listed", [
        ([("A", 35)], "A35H", "A35"),
        ([("A", 35), ("A", 72)], "A35H_A72H", "A35 A72"),
        ([("X", 16), ("X", 64)], "X16H_X64H", "X16 X64"),
    ])
    def test_tag_and_list(self, positions, tag, listed):
        assert his_tag(positions) == tag
        assert position_list(positions) == listed


class TestHisMeasurement:
    """One histidine, classified only when there is something to classify."""

    @pytest.mark.parametrize("pka,window,in_window", [
        (7.05, "ideal", True), (6.07, "usable", True), (5.35, "dead", False),
        (7.11, "ideal", True),
    ])
    def test_classifies_a_measured_pka(self, pka, window, in_window):
        site = HisMeasurement("A", 21, pka)
        assert site.window == window and site.in_window is in_window

    def test_an_unscored_histidine_has_no_window(self):
        """None is not 'dead' — PROPKA declining to answer is not an answer."""
        site = HisMeasurement("A", 72)
        assert site.window is None and site.in_window is False
        assert str(site) == "A72=NA"

    def test_renders_a_measurement_for_the_table(self):
        assert str(HisMeasurement("A", 21, 7.05)) == "A21=7.05"


class TestMeasureHis:
    """Mutate, write, measure — and keep the two failure modes apart."""

    def test_measures_every_requested_position(self, tmp_path):
        variant = measure_his(tmp_path / "parent.pdb", [("A", 21), ("A", 67)],
                              mutator=FakeMutator(), workdir=tmp_path,
                              propka=partial(fake_propka, table=MEASURED))
        assert [m.pka for m in variant.measurements] == [7.05, 7.11]
        assert variant.in_window == 2 and variant.ideal == 2
        assert variant.name == "parent__A21H_A67H"

    def test_a_failed_mutation_returns_none(self, tmp_path):
        """Distinct from a mutant PROPKA could not score — a broken build, not a hard question."""
        assert measure_his(tmp_path / "p.pdb", [("A", 21)],
                           mutator=FakeMutator(succeeds=False), workdir=tmp_path,
                           propka=partial(fake_propka, table=MEASURED)) is None

    def test_an_unscored_position_is_none_not_absent(self, tmp_path):
        """The silent failure: the variant exists, the histidine simply has no number."""
        variant = measure_his(tmp_path / "p.pdb", [("A", 21), ("A", 99)],
                              mutator=FakeMutator(), workdir=tmp_path,
                              propka=partial(fake_propka, table=MEASURED))
        assert len(variant.measurements) == 2
        assert variant.measurements[1].pka is None
        assert len(variant.measured) == 1
        assert variant.in_window == 1

    def test_mixed_windows_are_counted_separately(self, tmp_path):
        variant = measure_his(tmp_path / "p.pdb", [("A", 21), ("A", 35), ("A", 38)],
                              mutator=FakeMutator(), workdir=tmp_path,
                              propka=partial(fake_propka, table=MEASURED))
        assert variant.in_window == 2   # 7.05 ideal + 6.07 usable
        assert variant.ideal == 1       # only 7.05
        assert variant.measurements[2].window == "dead"

    def test_dest_keeps_the_mutant_outside_the_scratch_dir(self, tmp_path):
        """A combination that will be refolded has to outlive PROPKA's temp directory."""
        keep = tmp_path / "keep"
        keep.mkdir()
        mutator = FakeMutator()
        variant = measure_his(tmp_path / "p.pdb", [("A", 21)], mutator=mutator,
                              workdir=tmp_path, dest=keep,
                              propka=partial(fake_propka, table=MEASURED))
        assert variant.structure.parent == keep
        assert variant.structure.exists()

    def test_positions_round_trip(self, tmp_path):
        variant = measure_his(tmp_path / "p.pdb", [("A", 21), ("A", 67)],
                              mutator=FakeMutator(), workdir=tmp_path,
                              propka=partial(fake_propka, table=MEASURED))
        assert variant.positions == (("A", 21), ("A", 67))


class TestHisVariant:
    """Aggregates over one mutant's histidines."""

    def test_counts_are_zero_without_measurements(self):
        variant = HisVariant("parent", Path("p.pdb"),
                             (HisMeasurement("A", 21), HisMeasurement("A", 67)))
        assert variant.in_window == 0 and variant.ideal == 0 and variant.measured == ()


class TestHisCombinations:
    """Enumerating what to build, shared by the combination scan and the campaign builder."""

    def test_enumerates_every_subset_up_to_max_his(self):
        combos = his_combinations([("A", 21), ("A", 35), ("A", 67)], 3)
        assert len(combos) == 7  # 3 singles + 3 pairs + 1 triple
        assert [len(c) for c in combos] == [1, 1, 1, 2, 2, 2, 3]

    def test_max_his_below_the_site_count_truncates_by_size(self):
        combos = his_combinations([("A", 21), ("A", 35), ("A", 67)], 2)
        assert len(combos) == 6 and max(len(c) for c in combos) == 2

    def test_max_his_above_the_site_count_is_not_an_error(self):
        assert len(his_combinations([("A", 21)], 4)) == 1

    def test_cap_keeps_the_cheap_combinations(self):
        """Each added His costs interface confidence, so a cap must drop the large ones."""
        sites = [("A", 21), ("A", 35), ("A", 67), ("A", 70)]
        assert [len(c) for c in his_combinations(sites, 3, cap=5)] == [1, 1, 1, 1, 2]

    def test_no_sites_enumerates_nothing(self):
        assert his_combinations([], 3) == []

    def test_carries_whatever_shape_the_caller_uses(self):
        """The campaign builder keeps the measured pKa and verdict alongside each position."""
        scan_sites = [("A", 21, "R", 7.05, "ideal"), ("A", 35, "R", 6.07, "usable")]
        assert his_combinations(scan_sites, 2)[-1] == tuple(scan_sites)


class TestDeltaQ:
    """The charge a histidine gains between the two assay pH values."""

    def test_is_maximal_at_the_midpoint(self):
        """Why distance from 6.7 is the cost, not a low pKa."""
        assert delta_q(MIDPOINT) > max(delta_q(MIDPOINT - 0.4), delta_q(MIDPOINT + 0.4))

    @pytest.mark.parametrize("offset", [0.3, 0.9, 1.5])
    def test_penalises_both_directions_equally(self, offset):
        """A pKa above the midpoint is exactly as bad a bet as one the same distance below."""
        assert delta_q(MIDPOINT - offset) == pytest.approx(delta_q(MIDPOINT + offset))

    @pytest.mark.parametrize("pka,ph,fraction", [
        (6.7, 6.7, 0.5),      # half protonated at its own pKa
        (6.0, 8.0, 0.01),     # two units above: essentially neutral
        (7.4, 5.4, 0.99),     # two units below: essentially charged
    ])
    def test_henderson_hasselbalch(self, pka, ph, fraction):
        assert protonated_fraction(pka, ph) == pytest.approx(fraction, abs=0.005)


class TestSampledHis:
    """One histidine across the diffusion samples of its complex."""

    def test_floor_takes_the_sample_furthest_from_the_midpoint(self):
        """Not the lowest pKa: 7.55 is further from 6.7 than 6.61, so 7.55 sets the floor."""
        site = SampledHis("X", 21, (6.61, 7.55))
        assert site.delta_q_floor == pytest.approx(delta_q(7.55))
        assert site.delta_q_floor < delta_q(6.61)

    @pytest.mark.parametrize("pkas,in_window,robust", [
        ((6.61, 7.55), 2, True),          # both ends inside the 6.0-7.6 window
        ((3.99, 6.70), 1, False),         # helix0_R21H: structure-luck, the headline refutation
        ((5.56, 6.08), 1, False),         # helix0_E70H
        ((6.23, 6.50), 2, True),          # helix0_R72H
    ])
    def test_robustness_needs_every_sample(self, pkas, in_window, robust):
        site = SampledHis("X", 72, pkas)
        assert site.in_window_count == in_window and site.robust is robust

    def test_a_histidine_with_no_samples_scores_zero_rather_than_raising(self):
        """PROPKA declining to score every sample is a represented state: min() over nothing
        would take the whole ranking down instead of scoring that histidine 0."""
        empty = SampledHis("X", 21, ())
        assert empty.robust is False and empty.delta_q_floor == 0.0 and empty.spread == 0.0

    def test_spread_is_the_disagreement_between_samples(self):
        """2.71 units on one histidine of one design — the finding that forced the re-rank."""
        assert SampledHis("X", 21, (3.99, 6.70)).spread == pytest.approx(2.71)


class TestSwitchProfile:
    """Per-design aggregate — the number the submission was ranked on."""

    @pytest.mark.parametrize("design", list(PUBLISHED))
    def test_reproduces_the_submitted_delta_q_floor(self, design):
        profile = switch_profiles(robustness_rows(design))[0]
        assert profile.delta_q_floor == pytest.approx(PUBLISHED[design][1], abs=0.005)

    def test_a_parent_carrying_no_histidines_scores_zero(self):
        """tnf2_10 is the control and has to come out at exactly 0.00, not merely small."""
        assert SwitchProfile("free1_parent", ()).delta_q_floor == 0.0

    def test_counts_robust_histidines_as_the_table_reports_them(self):
        profile = switch_profiles(robustness_rows("helix0_R21H_E67H_E70H"))[0]
        assert profile.robust == 3 and profile.robust_tag == "3/3"

    def test_one_weak_histidine_does_not_hide_behind_two_strong_ones(self):
        """Summing lets a 2/3 design out-score a 1/1 -- how tnf2_02 beat tnf2_07."""
        weak = SwitchProfile("x", (SampledHis("X", 21, (4.64, 6.42)),
                                   SampledHis("X", 67, (7.27, 7.34)),
                                   SampledHis("X", 72, (6.36, 6.60))))
        strong = switch_profiles(robustness_rows("helix0_R72H"))[0]
        assert weak.robust == 2 and weak.delta_q_floor > strong.delta_q_floor

    def test_groups_rows_by_design_and_orders_histidines_by_position(self):
        rows = robustness_rows("helix0_R21H_E67H_E70H") + robustness_rows("helix0_R72H")
        profiles = {p.design: p for p in switch_profiles(rows)}
        assert len(profiles) == 2
        assert [h.resnum for h in profiles["helix0_R21H_E67H_E70H"].histidines] == [21, 67, 70]
