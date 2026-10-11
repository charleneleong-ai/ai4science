"""Reading PROPKA output, where the failure mode is silence rather than an error.

PROPKA omits a residue it could not score from the SUMMARY block without saying so. That is
how a whole scan once came back with every pKa missing: histidines had been renamed without
building the imidazole, so there was nothing to score. The parser must therefore be unable to
invent a value, and an absent residue must read as absent rather than as a default.
"""

from __future__ import annotations

from tnf_binder.pka import PKA_LINE, parse_summary, propka_executable

# A real SUMMARY block: two histidines scored, and ASP 11 alongside them.
SUMMARY = """\
propka3.5.1                                                        2026-10-11

---------  -----  ------   ---------------------------------------------------

SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP  11 A     3.80       3.80
   HIS  35 A     6.07       6.50
   HIS  72 A     6.66       6.50
   GLU  67 A     4.21       4.50
--------------------------------------------------------------------------------------------
Free energy of folding (kcal/mol) as a function of pH
"""

# The same prediction with HIS 72 missing, which is what a side-chain-less histidine looks like.
SUMMARY_MISSING_HIS = SUMMARY.replace("   HIS  72 A     6.66       6.50\n", "")


class TestParseSummary:
    """Only what PROPKA actually reported, keyed by chain and residue number."""

    def test_reads_every_scored_residue(self):
        pkas = parse_summary(SUMMARY)
        assert pkas[("A", 35)] == 6.07
        assert pkas[("A", 72)] == 6.66
        assert pkas[("A", 11)] == 3.80
        assert len(pkas) == 4

    def test_an_unscored_residue_is_absent_not_defaulted(self):
        """The silent failure: propka drops it, so the caller must see None, never a number."""
        pkas = parse_summary(SUMMARY_MISSING_HIS)
        assert ("A", 72) not in pkas
        assert pkas.get(("A", 72)) is None
        assert pkas[("A", 35)] == 6.07          # the rest still parse

    def test_stops_at_the_rule_after_the_table(self):
        """Text past the table must not be scraped for anything that looks like a pKa line."""
        trailing = SUMMARY + "   HIS 999 A     9.99       6.50\n"
        assert ("A", 999) not in parse_summary(trailing)

    def test_nothing_before_the_summary_header_is_read(self):
        assert parse_summary(SUMMARY.replace("SUMMARY OF THIS PREDICTION", "")) == {}

    def test_empty_output_is_empty(self):
        assert parse_summary("") == {}

    def test_negative_pka_parses(self):
        """Buried sites really do predict below zero; the regex has to accept the sign."""
        text = SUMMARY.replace("   HIS  35 A     6.07       6.50",
                               "   HIS  35 A    -0.05       6.50")
        assert parse_summary(text)[("A", 35)] == -0.05

    def test_distinguishes_chains_at_the_same_residue_number(self):
        """A trimer repeats residue numbers across protomers; collapsing them would be silent."""
        text = SUMMARY.replace("   HIS  72 A     6.66       6.50",
                               "   HIS  72 A     6.66       6.50\n   HIS  72 B     5.10       6.50")
        pkas = parse_summary(text)
        assert pkas[("A", 72)] == 6.66 and pkas[("B", 72)] == 5.10


class TestPkaLine:
    """The line pattern, which decides what counts as a scored residue."""

    def test_matches_an_ionizable_residue(self):
        assert PKA_LINE.match("   HIS  35 A     6.07       6.50")

    def test_ignores_a_non_ionizable_residue(self):
        assert PKA_LINE.match("   ALA  35 A     6.07       6.50") is None

    def test_ignores_the_header_row(self):
        assert PKA_LINE.match("       Group      pKa  model-pKa   ligand atom-type") is None


class TestPropkaExecutable:
    """Resolved next to the interpreter, because it is not on the subprocess PATH."""

    def test_returns_a_usable_name(self):
        assert propka_executable().endswith("propka3")
