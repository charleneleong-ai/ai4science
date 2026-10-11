"""The CSV hand-off between pipeline stages, and the coercion that has to fail loudly enough.

`number()` exists because BoltzGen and PROPKA both leave cells blank when a metric could not
be computed, and one blank must not sink a batch. But a silent default is the exact shape of
the worst bug in this project: a missing sequence identity defaulted to 0, which made the
fold-novelty gate look calibrated when it was not, and that claim had to be retracted. So the
default is tested explicitly, including that a caller can choose one that is not zero.
"""

from __future__ import annotations

import pytest

from tnf_binder.tables import design_stem, number, read_rows, write_rows


class TestReadRows:
    """A stage that has not run yet is the normal case, not an error."""

    def test_missing_file_is_an_empty_table(self, tmp_path):
        assert read_rows(tmp_path / "never_written.csv") == []

    def test_none_path_is_an_empty_table(self):
        assert read_rows(None) == []

    def test_header_only_file_is_an_empty_table(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("design,score\n")
        assert read_rows(path) == []

    def test_reads_rows_as_strings(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("design,score\nd1,0.5\n")
        assert read_rows(path) == [{"design": "d1", "score": "0.5"}]


class TestNumber:
    """Blank and malformed cells coerce; the default is the caller's choice, never a surprise."""

    def test_parses_a_normal_value(self):
        assert number({"x": "1.5"}, "x") == 1.5

    def test_blank_cell_takes_the_default(self):
        assert number({"x": ""}, "x") == 0.0

    def test_missing_key_takes_the_default(self):
        assert number({}, "x") == 0.0

    def test_junk_takes_the_default(self):
        assert number({"x": "n/a"}, "x") == 0.0

    def test_none_takes_the_default(self):
        assert number({"x": None}, "x") == 0.0

    def test_an_explicit_default_is_used_instead_of_zero(self):
        """A metric whose absence should not read as 'worst possible' needs a different floor."""
        for cell in ({"x": ""}, {}, {"x": "junk"}):
            assert number(cell, "x", default=-1.0) == -1.0

    def test_a_measured_zero_survives(self):
        """A genuine 0.0 must not read as absent — that is how a real measurement disappears."""
        assert number({"x": "0.0"}, "x", default=99.0) == 0.0
        assert number({"x": "0"}, "x", default=99.0) == 0.0

    def test_a_non_string_zero_hits_the_default(self):
        """The known sharp edge: `cell or default` swallows a numeric 0. read_rows only ever
        yields strings, so this bites a caller who builds rows itself."""
        assert number({"x": 0}, "x", default=99.0) == 99.0

    def test_negative_and_exponent_forms_parse(self):
        assert number({"x": "-0.5"}, "x") == -0.5
        assert number({"x": "1e-3"}, "x") == 0.001


class TestWriteRows:
    """Round-trips, and keeps the column order a reader will expect."""

    def test_round_trips(self, tmp_path):
        path = tmp_path / "t.csv"
        rows = [{"design": "d1", "score": "0.5"}, {"design": "d2", "score": "0.9"}]
        write_rows(path, rows)
        assert read_rows(path) == rows

    def test_explicit_columns_order_and_drop_extras(self, tmp_path):
        path = tmp_path / "t.csv"
        write_rows(path, [{"b": "2", "a": "1", "junk": "x"}], columns=["a", "b"])
        assert path.read_text().splitlines()[0] == "a,b"
        assert read_rows(path) == [{"a": "1", "b": "2"}]

    def test_no_rows_still_writes_a_readable_file(self, tmp_path):
        path = tmp_path / "t.csv"
        write_rows(path, [], columns=["a", "b"])
        assert path.read_text() == "a,b\r\n" or path.read_text() == "a,b\n"
        assert read_rows(path) == []


class TestDesignStem:
    """Reading the `design` column back to the design it names."""

    @pytest.mark.parametrize("name,expected", [
        ("F__helix0_R21H_R35H_E67H", "helix0_R21H_R35H_E67H"),
        ("FM__helix0_R72H", "helix0_R72H"),
        ("S2__helix0_E67H", "helix0_E67H"),
        ("P__screen_pocket_helix_0", "screen_pocket_helix_0"),
        ("helix0_R21H", "helix0_R21H"),           # already bare
    ])
    def test_strips_the_batch_marker(self, name, expected):
        assert design_stem(name) == expected

    def test_a_variant_keeps_its_his_tag(self):
        """Splitting on the LAST separator gave `A73H`, which names no design at all."""
        assert design_stem("V__screen_pocket_free_1__A73H") == "screen_pocket_free_1__A73H"

    def test_two_parents_sharing_a_his_tag_stay_distinct(self):
        """They collapsed onto one key, quietly dropping a row from the funnel."""
        assert (design_stem("V__screen_pocket_free_1__A73H")
                != design_stem("V__screen_pocket_free_2__A73H"))
