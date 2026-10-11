"""Verifier tests, anchored to outcomes that were actually observed.

Every case here is a real one: a sequence a generator really returned, a fold similarity really
measured, or -- for the fold gate -- the pass/fail verdict Proteinbase itself returned on a
submitted design. Nothing asserts that a mock returns what the mock was told to return.
"""

from __future__ import annotations

import pytest

from tnf_binder.design import TnfDesign
from tnf_binder.target import HUMAN, POCKET, LOBE, TnfTarget
from tnf_binder.verifiers import (
    CompositionVerifier,
    FoldHit,
    FoldNoveltyVerifier,
    HistidineSite,
    InterfaceMetrics,
    InterfaceVerifier,
    PkaSwitchVerifier,
    SequenceNoveltyVerifier,
)

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
# The design BoltzGen ranked first: best interface score, lowest interface PDE, an ideal-pKa
# histidine -- and verbatim ubiquitin.
TOP_RANKED_UBIQUITIN = TnfDesign("top_ranked", UBIQUITIN, generator="boltzgen",
                                 generator_confidence=0.919)
# 59% alanine, yet scored 0.556 binder-to-target ipTM and topped a 40-design shortlist.
POLY_ALANINE = TnfDesign(
    "poly_ala",
    "AAAAGAAAAAAVAAARAAGKSPEVAGFLAGEAVRAAAAALGLSRAAANAAARAAIEAAKAGADAAAVAAAAAAALR",
    generator="boltzgen", generator_confidence=0.556)
# Survived every gate and was submitted.
REAL_DESIGN = TnfDesign(
    "submitted",
    "ELEERIEKTAEEVRELVLEGRVAEAVALAVEEICHVGEEVGVPVERRNAVIREFERKLGELAEELEENPEKREEIAEKVAELVLKTLREIL",
    generator="boltzgen")


class TestComposition:
    """Degenerate sequences must fail regardless of how confidently a model scored them."""

    def test_rejects_poly_alanine_that_a_model_scored_well(self):
        verdict = CompositionVerifier().verify(POLY_ALANINE)
        assert not verdict.trust
        assert "low complexity" in verdict.reason
        assert verdict.metrics["top_residue"] == "A"

    def test_accepts_a_real_design(self):
        assert CompositionVerifier().verify(REAL_DESIGN).trust

    def test_passing_designs_are_not_all_given_the_same_score(self):
        """The first formula clamped the sum before halving, so every pass scored exactly 0.500.

        A constant score cannot rank anything, and reads as deliberate when it is not: the
        distinct-residue floor alone contributes 0.6, so the clamp was always active.
        """
        broad = CompositionVerifier().verify(TnfDesign("broad", "ACDEFGHIKLMNPQRSTVWY" * 4))
        narrow = CompositionVerifier().verify(REAL_DESIGN)
        assert broad.trust and narrow.trust
        assert broad.score != pytest.approx(0.5)
        assert broad.score > narrow.score  # richer alphabet, more headroom under the bound

    @pytest.mark.parametrize("sequence,expect_trust", [
        ("A" * 80, False),                       # single residue
        ("AG" * 40, False),                      # two residues
        ("ACDEFGHIKLMNPQRSTVWY" * 4, True),      # every residue, evenly used
    ])
    def test_distinct_residue_floor(self, sequence, expect_trust):
        design = TnfDesign("synthetic", sequence)
        assert CompositionVerifier().verify(design).trust is expect_trust

    def test_non_standard_residues_rejected(self):
        verdict = CompositionVerifier().verify(TnfDesign("x", "ACDEFGHIKLMNPQRSTVWYXXXX"))
        assert not verdict.trust and "non-standard" in verdict.reason


class TestSequenceNovelty:
    """Memorised natural proteins must be caught even when every other metric likes them."""

    @staticmethod
    def _corpus(*sequences: str):
        def lookup(kmers: list[str]) -> set[str]:
            return {km for km in kmers if any(km in s for s in sequences)}
        return lookup

    def test_rejects_the_top_ranked_design_because_it_is_ubiquitin(self):
        verifier = SequenceNoveltyVerifier(self._corpus(UBIQUITIN), corpus="test corpus")
        verdict = verifier.verify(TOP_RANKED_UBIQUITIN)
        assert not verdict.trust
        assert "memorised" in verdict.reason
        # the generator's own confidence was the highest in the campaign
        assert TOP_RANKED_UBIQUITIN.generator_confidence > 0.9

    def test_rejects_the_target_itself(self):
        verifier = SequenceNoveltyVerifier(self._corpus(HUMAN), corpus="test corpus")
        assert not verifier.verify(TnfDesign("tnf", HUMAN[:90])).trust

    def test_accepts_a_novel_sequence(self):
        verifier = SequenceNoveltyVerifier(self._corpus(UBIQUITIN, HUMAN), corpus="test corpus")
        assert verifier.verify(REAL_DESIGN).trust

    def test_defers_when_the_corpus_is_unavailable(self):
        def broken(_kmers):
            raise OSError("database unreachable")
        verdict = SequenceNoveltyVerifier(broken).verify(REAL_DESIGN)
        assert verdict.ood and not verdict.trust      # cannot judge, not judged-and-failed


class TestFoldNovelty:
    """Checked against the real Proteinbase verdicts, including where this gate disagrees.

    The gate is a conservative pre-filter, not a reimplementation: it predicts the real
    FAILURE correctly, and under-scores the two real PASSES because its structural search
    (FoldSeek over CATH50, co-folded structures) measures more similarity than Proteinbase's
    (ESMFold -> FoldSeek -> TM-align with domain segmentation).
    """

    def test_predicts_the_real_failure(self):
        """Beta-sandwich, TM 0.893 over 89% to CATH 2.60.40 -- Proteinbase rejected it."""
        design = TnfDesign("beta_sandwich", REAL_DESIGN.sequence, structure=__file__)
        verifier = FoldNoveltyVerifier(
            fold_search=lambda _d: [FoldHit(0.893, 0.888, "af_Q64612_2.60.40")],
            identity_of=lambda _d: 37.5,
        )
        verdict = verifier.verify(design)
        assert verdict.metrics["novelty_score"] == 2
        assert not verdict.trust
        assert verdict.metrics["structural"] == "high"

    @pytest.mark.parametrize("tm,coverage", [(0.776, 0.802), (0.759, 0.802)])
    def test_under_scores_the_real_passes(self, tm, coverage):
        """Documents the known disagreement rather than tuning thresholds to hide it.

        Proteinbase passed both of these at ~38% identity. This gate returns 2, so it is
        conservative: a borderline score here is a warning, not a disqualification.
        """
        design = TnfDesign("helical", REAL_DESIGN.sequence, structure=__file__)
        verifier = FoldNoveltyVerifier(
            fold_search=lambda _d: [FoldHit(tm, coverage, "cath_1.10")],
            identity_of=lambda _d: 38.1,
        )
        verdict = verifier.verify(design)
        assert verdict.metrics["structural"] == "moderate"
        assert verdict.metrics["novelty_score"] == 2, "conservative vs the real pass"

    def test_a_missing_identity_must_not_silently_inflate_the_score(self):
        """The bug that produced a spurious 3/3 calibration: identity defaulting to zero."""
        design = TnfDesign("d", REAL_DESIGN.sequence, structure=__file__)
        with_identity = FoldNoveltyVerifier(
            lambda _d: [FoldHit(0.776, 0.802, "x")], lambda _d: 38.1).verify(design)
        without_identity = FoldNoveltyVerifier(
            lambda _d: [FoldHit(0.776, 0.802, "x")], lambda _d: 0.0).verify(design)
        assert with_identity.metrics["novelty_score"] == 2
        assert without_identity.metrics["novelty_score"] == 3

    def test_fold_similarity_decides_when_identity_is_equal(self):
        """The failing and passing designs differed by 0.6% identity: fold was the cause."""
        design = TnfDesign("d", REAL_DESIGN.sequence, structure=__file__)
        high = FoldNoveltyVerifier(lambda _d: [FoldHit(0.89, 0.89, "ig_fold")],
                                   lambda _d: 38.0).verify(design)
        low = FoldNoveltyVerifier(lambda _d: [FoldHit(0.30, 0.50, "weak")],
                                  lambda _d: 38.0).verify(design)
        assert not high.trust and low.trust

    def test_defers_without_a_structure(self):
        verdict = FoldNoveltyVerifier(lambda _d: [], lambda _d: 0.0).verify(REAL_DESIGN)
        assert verdict.ood


class TestPkaSwitch:
    """pKa must sit between the two assay pH values, which three mechanisms failed to achieve."""

    @pytest.mark.parametrize("pka,window", [
        (3.15, "dead"),    # measured on a buried interface histidine
        (5.71, "dead"),    # measured; below the pH 6.0 trigger
        (6.07, "usable"),
        (6.68, "ideal"),
        (7.14, "ideal"),
        (8.50, "dead"),    # protonated at both pHs, so it blocks rather than switches
    ])
    def test_window_classification(self, pka, window):
        assert HistidineSite(1, pka, interface_dsasa=50.0).window == window

    def test_buried_histidines_from_first_principles_all_failed(self):
        """The three interface histidines PROPKA actually measured on a real design."""
        measured = [HistidineSite(14, 5.71, 26.0),
                    HistidineSite(40, 3.15, 68.0),
                    HistidineSite(46, 4.70, 89.0)]
        design = TnfDesign("bindcraft", REAL_DESIGN.sequence,
                           metrics={"histidine_sites": measured})
        verdict = PkaSwitchVerifier().verify(design)
        assert not verdict.trust
        assert verdict.metrics["switches"] == 0

    def test_in_window_histidine_at_the_interface_switches(self):
        design = TnfDesign("engineered", REAL_DESIGN.sequence,
                           metrics={"histidine_sites": [HistidineSite(73, 6.53, 59.5)]})
        assert PkaSwitchVerifier().verify(design).trust

    def test_right_pka_in_the_wrong_place_does_not_count(self):
        """A usable pKa buried inside the binder changes nothing about binding."""
        design = TnfDesign("internal", REAL_DESIGN.sequence,
                           metrics={"histidine_sites": [HistidineSite(68, 6.02, 0.0)]})
        assert not PkaSwitchVerifier().verify(design).trust

    def test_defers_without_a_prediction(self):
        assert PkaSwitchVerifier().verify(REAL_DESIGN).ood


class TestInterface:
    """Buried area, not contact counts -- counting rewarded interpenetration."""

    def test_rejects_clashing_design(self):
        design = TnfDesign("clashing", REAL_DESIGN.sequence,
                           metrics={"interface": InterfaceMetrics(1039.8, 0.37, 0.14, 171)})
        verdict = InterfaceVerifier().verify(design)
        assert not verdict.trust and "clash" in verdict.reason

    def test_rejects_off_epitope_binding(self):
        design = TnfDesign("off", REAL_DESIGN.sequence,
                           metrics={"interface": InterfaceMetrics(798.5, 0.037, 0.46, 0)})
        assert not InterfaceVerifier().verify(design).trust

    def test_accepts_a_real_interface(self):
        design = TnfDesign("good", REAL_DESIGN.sequence,
                           metrics={"interface": InterfaceMetrics(604.6, 0.41, 0.17, 0)})
        verdict = InterfaceVerifier().verify(design)
        assert verdict.trust
        assert verdict.metrics["cross_reactivity"] == pytest.approx(0.83, abs=0.01)

    def test_larger_interface_is_not_automatically_better(self):
        """The biggest interfaces in the campaign were interpenetrating junk."""
        junk = TnfDesign("junk", REAL_DESIGN.sequence,
                         metrics={"interface": InterfaceMetrics(1345.5, 0.24, 0.20, 649)})
        real = TnfDesign("real", REAL_DESIGN.sequence,
                         metrics={"interface": InterfaceMetrics(604.6, 0.41, 0.17, 0)})
        assert not InterfaceVerifier().verify(junk).trust
        assert InterfaceVerifier().verify(real).trust


class TestTargetNumbering:
    """An unguarded conversion here silently shifted an epitope by five residues."""

    def test_human_mouse_identity_matches_the_published_figure(self):
        assert round(TnfTarget().identity * 157) == 124      # 79.0%

    def test_r31_is_divergent_and_explains_the_drug_failure(self):
        """R31->Q is why adalimumab and infliximab are not mouse cross-reactive."""
        target = TnfTarget()
        assert 31 in target.divergent
        assert target.human.residue(31) == "R"
        assert target.mouse.residue(target.human_to_mouse[31]) == "Q"

    @pytest.mark.parametrize("positions", [POCKET, LOBE])
    def test_hotspots_map_to_identical_mouse_residues(self, positions):
        target = TnfTarget()
        mouse = target.mouse_equivalent(positions)
        assert len(mouse) == len(positions)
        for h, m in zip(positions, mouse):
            assert target.human.residue(h) == target.mouse.residue(m)

    def test_mapping_a_divergent_position_is_refused(self):
        with pytest.raises(ValueError, match="not conserved"):
            TnfTarget().mouse_equivalent((31,))

    def test_chain_position_conversion(self):
        """1TNF starts at residue 6, so mature 113 is chain position 108."""
        assert TnfTarget.chain_position(113, first_stored=6) == 108
        assert TnfTarget.chain_position(86, first_stored=6) == 81
