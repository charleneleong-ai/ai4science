"""The admissibility gate, tested against copies it has to catch.

The references are the real chains of 7OW6 at their deposited lengths, and the design is
bulge_core_23 from the first BoltzGen sweep. Both of those matter. The gate's failure mode was
that a truncated copy of a target chain scored zero, and how badly the superseded measure
underestimates depends on the ratio of the lengths: a 40-residue slice of the 276-residue heavy
chain scored 0.000, while the same slice against a 120-residue excerpt scores 0.47 and looks
caught. A shortened reference hides the bug the test exists to pin down.
"""

from __future__ import annotations

import pytest
from novelty import (
    DEFAULT_THRESHOLD,
    copied_fraction,
    global_aligner,
    global_identity,
    judge,
    local_aligner,
)

HLA = (  # 7OW6 chain A, HLA-A*11:01 heavy chain
    "GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQS"
    "EDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGK"
    "ETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHE"
    "GLPKPLTLRWEP"
)
B2M = (  # 7OW6 chain B, the chain the generator actually hands back
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLS"
    "QPKIVKWDRDM"
)
TCR_BETA = (  # 7OW6 chain E, the known receptor, which the rules forbid modifying
    "MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHEYMYWYRQDPGMGLRLIHYSVGEGTTAKGEVPDGYNVSRLKKQNFLLGLESAAPSQTS"
    "VYFCASKVGPGQHNSPLHFGNGTRLTVTEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDP"
    "QPLKEQPALNDSRYALSSRLRVSATFWQDPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRA"
)
EPITOPE_G12D = "VVVGADGVGK"
DESIGN = "EAAREELFRLALELLRKTDADEAFERVVERAKELGLSEEEAAAIVGEASRVLLREREERER"

REFS = {
    "hla_heavy": HLA,
    "b2m": B2M,
    "tcr_E": TCR_BETA,
    "peptide_g12d": EPITOPE_G12D,
}


@pytest.fixture(scope="module")
def local():
    return local_aligner()


def substitute(sequence: str, every: int) -> str:
    """Replace one residue in `every`, with one that is not already there."""
    out = list(sequence)
    for i in range(0, len(out), every):
        out[i] = "W" if out[i] != "W" else "P"
    return "".join(out)


class TestCopiedFraction:
    """A design that is part of a target chain has to read as copied."""

    @pytest.mark.parametrize(
        "name,sequence,floor",
        [
            ("b2m_verbatim", B2M, 1.0),
            ("b2m_substituted", substitute(B2M, 10), 0.85),
            ("hla_fragment_40", HLA[50:90], 1.0),
            ("hla_fragment_80", HLA[50:130], 1.0),
            ("hla_fragment_substituted", substitute(HLA[50:130], 20), 0.90),
            ("tcr_fragment_60", TCR_BETA[40:100], 1.0),
        ],
    )
    def test_copies_are_caught(self, name, sequence, floor, local):
        worst = max(copied_fraction(sequence, ref, local) for ref in REFS.values())
        assert worst >= floor, f"{name} scored {worst:.3f}, under {floor}"
        assert worst >= DEFAULT_THRESHOLD

    def test_de_novo_design_passes(self, local):
        worst = max(copied_fraction(DESIGN, ref, local) for ref in REFS.values())
        assert worst < DEFAULT_THRESHOLD, f"real design scored {worst:.3f}"

    def test_normalised_by_the_design_not_the_reference(self, local):
        """Against the 10-mer epitope, chance matches must not become a large fraction."""
        assert copied_fraction(DESIGN, EPITOPE_G12D, local) < 0.1

    def test_empty_sequence_is_not_a_copy(self, local):
        assert copied_fraction("", HLA, local) == 0.0
        assert copied_fraction(HLA, "", local) == 0.0


class TestSupersededGlobalMetric:
    """Locked in because this is the hole the gate shipped with, and it is invisible to review."""

    @pytest.mark.parametrize(
        "fragment,reference", [(HLA[50:90], HLA), (TCR_BETA[40:100], TCR_BETA)]
    )
    def test_global_identity_cannot_see_a_truncated_copy(self, fragment, reference):
        assert global_identity(fragment, reference, global_aligner()) < DEFAULT_THRESHOLD


class TestJudge:
    def test_rejects_the_copy_and_keeps_the_design(self):
        verdicts = {
            v.design: v
            for v in judge(
                {"copy": HLA[50:130], "design": DESIGN}, REFS, DEFAULT_THRESHOLD
            )
        }
        assert not verdicts["copy"].admissible
        assert verdicts["copy"].worst_match == "hla_heavy"
        assert verdicts["design"].admissible

    def test_whole_chain_identity_is_still_reported(self):
        """The superseded number stays in the row, so a reader can see why the metric changed."""
        [verdict] = judge({"copy": HLA[50:130]}, REFS, DEFAULT_THRESHOLD)
        assert verdict.row["copied_fraction"] == 1.0
        assert verdict.row["whole_chain_identity"] < DEFAULT_THRESHOLD
