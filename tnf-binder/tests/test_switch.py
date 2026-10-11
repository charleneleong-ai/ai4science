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
    HisMeasurement,
    HisVariant,
    his_tag,
    measure_his,
    position_list,
)

# measured on screen_pocket_helix_0: two ideal, one usable, one dead
MEASURED = {("A", 21): 7.05, ("A", 67): 7.11, ("A", 35): 6.07, ("A", 38): 5.35}


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
