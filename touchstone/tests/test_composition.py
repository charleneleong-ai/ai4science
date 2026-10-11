from pathlib import Path

import numpy as np
import pytest

from touchstone import BinderDesign, CompositionVerifier, rank
from touchstone.core import CoordinationSite
from touchstone.service import verify_structure

PACKED = Path(__file__).parent / "fixtures" / "ligmpnn_nickel_packed.pdb"

# human ubiquitin (1-76) — the positive control: a real protein must clear the gate
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
# the campaign's worst offender rebuilt: 59% alanine (0.556 binder-to-target ipTM, top of a
# 40-design shortlist) but 20 distinct residues, so only the single-residue bound can catch it
ALANINE_RICH = "A" * 27 + "GSTLKEWVNDPQRCYFMIH"


def _design(seq: str) -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
    return BinderDesign(seq, site, generator="test", generator_confidence=0.5)


@pytest.fixture
def verifier() -> CompositionVerifier:
    return CompositionVerifier()


class TestDegenerateSequences:
    """The known negatives from the 57-design campaign: low-complexity Ala/Gly runs."""

    @pytest.mark.parametrize(
        "seq, expected_reason",
        [
            ("A" * 60, "low complexity"),
            ("G" * 60, "low complexity"),
            (ALANINE_RICH, "low complexity"),
            ("G" * 22 + "CDEFHIKLMNPQRSTVWY" * 2, "low complexity"),  # the 37%-glycine design
            ("AAAAAGGGGG" + "CDEFHIKLMN", "Ala+Gly"),  # 25% A + 25% G, under every other bound
            ("CDEFHIKLMNP" * 4, "distinct residues"),  # 11 distinct, no residue over 10%
        ],
        ids=["poly_ala", "poly_gly", "ala_59pct", "gly_37pct", "ala_gly_half", "low_diversity"],
    )
    def test_rejected_as_weak(self, verifier, seq, expected_reason):
        v = verifier.verify(_design(seq))
        assert v.label == "weak" and v.score == 0.0
        assert expected_reason in v.reason

    def test_metrics_expose_the_offending_residue(self, verifier):
        m = verifier.verify(_design(ALANINE_RICH)).metrics
        assert m["top_residue"] == "A" and m["top_fraction"] == pytest.approx(0.587, abs=0.005)
        assert m["distinct"] == 20  # the distinct-count bound can't catch it — composition is the catch


class TestJudgeability:
    """defer means 'cannot judge'; weak means 'judged and failed'. Load-bearing distinction."""

    def test_empty_sequence_defers(self, verifier):
        v = verifier.verify(_design(""))
        assert v.label == "defer" and "empty sequence" in v.reason

    def test_nonstandard_residues_are_weak_not_defer(self, verifier):
        v = verifier.verify(_design("X" + UBIQUITIN[1:]))
        assert v.label == "weak" and "non-standard residues ['X']" in v.reason


class TestRealSequences:
    def test_real_protein_trusts(self, verifier):
        assert verifier.verify(_design(UBIQUITIN)).label == "trust"

    def test_score_prefers_diversity_on_the_shared_zero_to_one_scale(self, verifier):
        borderline = "K" * 7 + "CDEFHILMNPQRSTVWY"  # 29% K, just inside max_single_fraction
        good, meh = (verifier.verify(_design(s)).score for s in (UBIQUITIN, borderline))
        assert 0.0 < meh < good <= 1.0  # no saturation: it must grade, not emit a constant

    def test_thresholds_are_tunable(self):
        strict = CompositionVerifier(max_single_fraction=0.10)
        assert not strict.verify(_design(UBIQUITIN)).trust  # 11.8% Leu now trips the bound


def test_ranks_a_degenerate_design_below_a_real_one():
    # the opt-in wiring: a Verifier the pipeline can rank with, judging sequence alone
    ranked = rank([_design(ALANINE_RICH), _design(UBIQUITIN)], CompositionVerifier())
    assert [d.sequence for d, _ in ranked] == [UBIQUITIN, ALANINE_RICH]


class TestServiceWiring:
    """Opt-in, but discoverable: the tier is advertised in the stack and flipped by one flag."""

    def test_advertised_but_not_run_by_default(self):
        r = verify_structure(PACKED, "Ni2+")
        by = {s["stage"]: s["status"] for s in r["stack"]}
        assert by["composition"] == "needs_input"  # never silently folded into the default consensus
        assert "composition" in r["not_run"]

    def test_needs_a_sequence_not_just_the_flag(self):
        # the gate's only input is the sequence — without one it stays unrun rather than deferring
        r = verify_structure(PACKED, "Ni2+", composition=True)
        assert {s["stage"]: s["status"] for s in r["stack"]}["composition"] == "needs_input"

    @pytest.mark.parametrize("seq, label", [(UBIQUITIN, "trust"), (ALANINE_RICH, "weak")])
    def test_flag_plus_sequence_runs_the_tier(self, seq, label):
        r = verify_structure(PACKED, "Ni2+", composition=True, sequence=seq)
        assert {s["stage"]: s["status"] for s in r["stack"]}["composition"] == "ran"
        assert r["verifiers"]["composition"]["label"] == label
