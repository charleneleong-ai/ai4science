"""Calibrating a gate against designs whose outcome is known.

The counts here are the real ones measured over 4065 Proteinbase designs carrying an
experimental `expressed` verdict, so a change that breaks the arithmetic shows up as a verdict
that contradicts something already observed rather than as a failed synthetic assertion.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from touchstone.calibration import (
    BINDS,
    EXPRESSED,
    Calibration,
    Outcome,
    bands,
    calibrate,
    load_export,
)

# Measured over the export: expression by cysteine count, within the 60-100 aa band.
ONE_CYS = (133, 127)      # n, expressed -> 95.5%
NO_CYS = (1018, 979)      # n, expressed -> 96.2%
# Odd and even non-zero cysteine counts, pooled over all lengths: identical rates.
ODD, EVEN = (646, 560), (1120, 971)


def population(groups: dict[bool, tuple[int, int]]) -> list[Outcome]:
    """Designs whose only distinguishing feature is whether a gate would keep them.

    The sequence encodes the gate's answer so a predicate can read it back without the test
    needing a real feature; the arithmetic under test is the same either way.
    """
    out = []
    for keep, (n, successes) in groups.items():
        tag = "K" if keep else "D"
        out += [Outcome(sequence=tag, expressed=i < successes) for i in range(n)]
    return out


def keeps_tagged(design: Outcome) -> bool:
    return design.sequence == "K"


class TestCalibration:
    """The arithmetic a gate's track record rests on."""

    def test_reports_both_rates_and_the_base(self):
        c = calibrate(population({True: NO_CYS, False: ONE_CYS}), keeps_tagged)
        assert c.kept == 1018 and c.dropped == 133
        assert c.kept_rate == pytest.approx(0.962, abs=0.001)
        assert c.dropped_rate == pytest.approx(0.955, abs=0.001)
        assert c.base_rate == pytest.approx(0.961, abs=0.001)

    def test_the_cysteine_gate_does_not_discriminate(self):
        """0.7 points over 1151 designs: the liability argued from first principles."""
        c = calibrate(population({True: NO_CYS, False: ONE_CYS}), keeps_tagged)
        assert c.separation == pytest.approx(0.007, abs=0.002)
        assert not c.discriminates and c.verdict == "undecided"

    def test_odd_versus_even_cysteine_is_exactly_nothing(self):
        """86.7% against 86.7%: the unpairable-cysteine argument, refuted on 1766 designs."""
        c = calibrate(population({True: EVEN, False: ODD}), keeps_tagged)
        assert c.kept_rate == pytest.approx(c.dropped_rate, abs=0.002)
        assert not c.discriminates

    def test_a_real_effect_does_discriminate(self):
        """Length: 96.8% at 60-80 aa against 79.6% above 150."""
        c = calibrate(population({True: (349, 338), False: (543, 432)}), keeps_tagged)
        assert c.separation > 0.15 and c.discriminates and c.verdict == "discriminates"

    def test_a_gate_that_keeps_everything_is_vacuous_not_perfect(self):
        """The expression verifier passed 13 of 13 at its shipped defaults and read as a pass."""
        c = calibrate(population({True: (13, 13)}), keeps_tagged)
        assert c.kept_rate == 1.0 and not c.discriminates and c.verdict == "vacuous"

    def test_a_backwards_gate_is_reported_as_inverted_not_as_failure(self):
        c = calibrate(population({True: (543, 432), False: (349, 338)}), keeps_tagged)
        assert c.separation < 0 and c.verdict == "inverted"

    def test_unmeasured_designs_are_excluded_not_counted_as_failures(self):
        measured = [Outcome("K", expressed=True), Outcome("D", expressed=False)]
        unmeasured = [Outcome("K"), Outcome("D")]
        assert calibrate(measured, keeps_tagged) == calibrate(measured + unmeasured, keeps_tagged)

    def test_an_outcome_nothing_measured_leaves_it_vacuous(self):
        """`binds` is absent from most of the export; asking for it must not read as 0%."""
        c = calibrate([Outcome("K", expressed=True)], keeps_tagged, outcome=BINDS)
        assert c.total == 0 and c.verdict == "vacuous"

    def test_separation_is_not_lift(self):
        """A gate keeping most of the population barely beats the base rate it is part of."""
        c = calibrate(population({True: (900, 880), False: (100, 60)}), keeps_tagged)
        assert c.separation > 0.3 and c.lift < 0.05


class TestBands:
    """Looking at a feature before choosing a threshold on it."""

    def test_rates_per_slice(self):
        designs = ([Outcome("x" * 70, expressed=i < 338) for i in range(349)]
                   + [Outcome("x" * 200, expressed=i < 432) for i in range(543)])
        got = bands(designs, lambda d: len(d.sequence), (60, 100, 250))
        assert [b.n for b in got] == [349, 543]
        assert got[0].rate == pytest.approx(0.968, abs=0.002)
        assert got[1].rate == pytest.approx(0.796, abs=0.002)

    def test_an_empty_slice_reports_zero_n_rather_than_a_rate(self):
        got = bands([Outcome("xx", expressed=True)], lambda d: len(d.sequence), (0, 1, 10))
        assert got[0].n == 0 and got[0].rate == 0.0 and got[1].n == 1


class TestLoadExport:
    """Reading Proteinbase's CSV without mistaking a prediction for a measurement."""

    @staticmethod
    def write(tmp_path: Path, evaluations: list[dict], sequence: str = "MKLV") -> Path:
        path = tmp_path / "export.csv"
        with path.open("w", newline="") as fh:
            fh.write('"id","name","sequence","author","designMethod","evaluations"\n')
            fh.write(f'"a","n","{sequence}","x","boltzgen","{json.dumps(evaluations).replace(chr(34), chr(34) * 2)}"\n')
        return path

    def test_reads_an_experimental_verdict(self, tmp_path):
        path = self.write(tmp_path, [{"type": "experimental", "metric": "expressed",
                                      "value": "true", "target": "tnf-a"}])
        assert load_export(path)[0].expressed is True

    def test_a_computational_expressed_row_is_not_a_measurement(self, tmp_path):
        """Calibrating a predictor against another predictor measures agreement, not truth."""
        path = self.write(tmp_path, [{"type": "computational", "metric": "expressed",
                                      "value": "true"}])
        assert load_export(path)[0].expressed is None

    def test_any_replicate_expressing_counts_as_expressed(self, tmp_path):
        path = self.write(tmp_path, [
            {"type": "experimental", "metric": "expressed", "value": "false"},
            {"type": "experimental", "metric": "expressed", "value": "true"}])
        assert load_export(path)[0].expressed is True

    def test_kd_takes_the_tightest_measurement(self, tmp_path):
        path = self.write(tmp_path, [{"type": "experimental", "metric": "kd", "value": "23.7"},
                                     {"type": "experimental", "metric": "kd", "value": "491"}])
        assert load_export(path)[0].kd_nm == pytest.approx(23.7)

    def test_target_filter_keeps_only_designs_measured_against_it(self, tmp_path):
        path = self.write(tmp_path, [{"type": "experimental", "metric": "expressed",
                                      "value": "true", "target": "egfr"}])
        assert load_export(path, target="egfr") and not load_export(path, target="tnf-a")

    def test_a_row_without_a_sequence_is_skipped(self, tmp_path):
        path = self.write(tmp_path, [{"type": "experimental", "metric": "expressed",
                                      "value": "true"}], sequence="")
        assert load_export(path) == []


class TestCalibrationRendering:
    def test_str_carries_both_groups_and_the_verdict(self):
        c = Calibration(kept=100, dropped=100, kept_successes=90, dropped_successes=50)
        rendered = str(c)
        assert "discriminates" in rendered and "90.0%" in rendered and "50.0%" in rendered

    def test_outcome_defaults_to_expression(self):
        assert Calibration(1, 1, 1, 0).outcome == EXPRESSED


class TestCalibrateCommand:
    """`touchstone calibrate` — the thing to run on day one of a new target."""

    @staticmethod
    def export(tmp_path: Path) -> Path:
        """Two groups that differ only in length, so the length gate must separate them."""
        path = tmp_path / "export.csv"
        rows = [("M" * 70, "true")] * 40 + [("M" * 200, "false")] * 40
        with path.open("w", newline="") as fh:
            fh.write('"id","name","sequence","author","designMethod","evaluations"\n')
            for i, (seq, call) in enumerate(rows):
                ev = f'[{{""type"":""experimental"",""metric"":""expressed"",""value"":""{call}""}}]'
                fh.write(f'"{i}","n","{seq}","x","m","{ev}"\n')
        return path

    def test_reports_a_verdict_per_gate(self, tmp_path):
        from typer.testing import CliRunner

        from touchstone.cli import app
        result = CliRunner().invoke(app, ["calibrate", str(self.export(tmp_path))])
        assert result.exit_code == 0
        assert "80 designs with a measured expressed" in result.stdout
        # Asserts on the lines outside the table: rich reflows and elides the table body at
        # the runner's terminal width, so its cells are not reliably capturable. The verdict
        # logic those cells render is unit-tested directly in TestCalibration.
        assert "60 to 80: 100.0% (n=40)" in result.stdout
        assert "150 to 10000: 0.0% (n=40)" in result.stdout

    def test_exits_nonzero_when_nothing_was_measured(self, tmp_path):
        from typer.testing import CliRunner

        from touchstone.cli import app
        empty = tmp_path / "empty.csv"
        empty.write_text('"id","name","sequence","author","designMethod","evaluations"\n'
                         '"a","n","MKLV","x","m","[]"\n')
        result = CliRunner().invoke(app, ["calibrate", str(empty)])
        assert result.exit_code == 1
