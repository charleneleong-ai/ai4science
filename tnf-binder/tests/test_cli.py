"""CLI behaviour: the submission pre-flight rules, and that each command renders and exits right.

The submission rules are Proteinbase's stated format, and `tnf_submission_v1.csv` is the file
actually uploaded -- it must keep passing.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from tnf_binder.cli import app, csv_problems, row_problem

SUBMISSION = Path(__file__).parents[1] / "submission" / "tnf_submission_v1.csv"
runner = CliRunner()


def row(name: str = "d1", sequence: str = "ACDEFGHIKLMNPQ",
        molecule_class: str = "single_chain") -> dict[str, str]:
    return {"name": name, "sequence": sequence, "molecule_class": molecule_class}


class TestRowProblem:
    """One submission row against the stated format."""

    def test_accepts_a_well_formed_row(self):
        assert row_problem(row()) is None

    @pytest.mark.parametrize("field,value,expected", [
        ("molecule_class", "minibinder", "bad class"),
        ("sequence", "ACDEFGHIKLMNPQX", "non-standard residues"),
        ("sequence", "ACDEF", "length 5"),
        ("sequence", "A" * 251, "length 251"),
    ])
    def test_rejects(self, field, value, expected):
        assert expected in row_problem(row(**{field: value}))

    def test_reports_class_before_sequence(self):
        """Class is checked first so an unsupported class is never masked by a length complaint."""
        assert "bad class" in row_problem(row(sequence="ACDEF", molecule_class="minibinder"))


class TestCsvProblems:
    """A whole submission file."""

    def test_the_uploaded_submission_still_passes(self):
        with SUBMISSION.open() as handle:
            assert csv_problems(list(csv.DictReader(handle))) == []

    def test_empty_file(self):
        assert csv_problems([]) == ["empty"]

    def test_duplicate_sequences(self):
        assert "duplicate sequences" in csv_problems([row("a"), row("b")])

    def test_row_cap(self):
        rows = [row(f"d{i}", "ACDEFGHIKLMNPQ" + "G" * i) for i in range(21)]
        assert "21 rows > 20 limit" in csv_problems(rows)

    def test_wrong_columns(self):
        problems = csv_problems([{"sequence": "ACDEFGHIKLMNPQ", "name": "d1",
                                  "molecule_class": "single_chain"}])
        assert any("columns" in p for p in problems)

    def test_collects_every_bad_row_not_just_the_first(self):
        rows = [row("a", "ACDEF"), row("b", "ACDEFGHIKLMNPQX")]
        assert len(csv_problems(rows)) == 2


class TestCommands:
    """Each command renders and exits with the right code."""

    def test_epitope_json_reports_the_computed_epitope(self):
        result = runner.invoke(app, ["epitope", "--json"])
        assert result.exit_code == 0
        payload = json.loads(result.stdout)
        assert payload["epitope_positions"] == 33
        assert len(payload["conserved"]) == 24
        assert payload["identity"] == pytest.approx(0.790, abs=1e-3)
        assert 31 in payload["divergent"]

    def test_epitope_table_flags_the_divergent_positions(self):
        """Position 31 is human R / mouse Q, and has to render as not conserved."""
        result = runner.invoke(app, ["epitope"])
        assert result.exit_code == 0
        row_31 = next(line for line in result.stdout.splitlines() if line.startswith("│ 31 "))
        assert [cell.strip() for cell in row_31.strip("│").split("│")] == ["31", "R", "Q", "NO"]
        assert "R31->Q is the trap" in result.stdout

    def test_validate_passes_the_real_submission_directory(self):
        result = runner.invoke(app, ["validate", str(SUBMISSION.parent)])
        assert result.exit_code == 0
        assert "OK" in result.stdout

    def test_validate_exits_nonzero_on_a_bad_file(self, tmp_path):
        (tmp_path / "bad.csv").write_text("name,sequence,molecule_class\nd1,ACDEF,single_chain\n")
        result = runner.invoke(app, ["validate", str(tmp_path)])
        assert result.exit_code == 1

    def test_confidence_exits_nonzero_when_nothing_is_folded(self, tmp_path):
        result = runner.invoke(app, ["confidence", str(tmp_path)])
        assert result.exit_code == 2

    def test_help_lists_every_subcommand(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        for command in ("epitope", "verify", "validate", "confidence", "wandb-log"):
            assert command in result.stdout

    def test_bare_invocation_demands_a_subcommand(self):
        """The root callback exists so running the bare app is an error, not a no-op."""
        assert runner.invoke(app, []).exit_code == 2
