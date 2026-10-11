"""Generator-input numbering, which has silently aimed a campaign at the wrong surface once.

BoltzGen's `parse_range` counts chain positions as stored in the file, not residue numbers.
1TNF starts at residue 6, so writing mature numbers straight into a spec shifts the epitope by
five. These tests build structures that start at an arbitrary residue number and assert the
guard both maps correctly and refuses when the residues are not the ones expected.
"""

from __future__ import annotations

from pathlib import Path

import itertools

import pytest
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.PDB import PDBIO, StructureBuilder

from tnf_binder.specs import (
    StoredNumbering,
    chain_block,
    cofold_spec,
    collapse_ranges,
    helix_bundle,
    secondary_structure,
)

FIRST_RESIDUE = 6  # as 1TNF stores it, the offset that caused the original shift


def build(path: Path, sequence: str, *, start: int = FIRST_RESIDUE, chain: str = "A") -> Path:
    """One chain whose first modelled residue is numbered `start`."""
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure("s")
    builder.init_model(0)
    builder.init_chain(chain)
    builder.init_seg(" ")
    for offset, residue in enumerate(sequence):
        builder.init_residue(protein_letters_1to3[residue].upper(), " ", start + offset, " ")
        builder.init_atom("CA", (float(offset), 0.0, 0.0), 0.0, 1.0, " ", "CA", element="C")
    io = PDBIO()
    io.set_structure(builder.get_structure())
    io.save(str(path))
    return path


@pytest.fixture
def numbering(tmp_path) -> StoredNumbering:
    # residue 6 is 'A', 7 'C', 8 'D' ... so mature n sits at chain position n - 5
    return StoredNumbering.from_file(build(tmp_path / "s.pdb", "ACDEFGHIKL"))


class TestCollapseRanges:
    """BoltzGen's 'a..b,c' syntax."""

    @pytest.mark.parametrize("numbers,expected", [
        ([1, 2, 3], "1..3"),
        ([1], "1"),
        ([1, 3, 5], "1,3,5"),
        ([1, 2, 3, 7, 8], "1..3,7..8"),
        ([3, 1, 2], "1..3"),          # unsorted input
        ([5, 5, 6], "5..6"),          # duplicates collapse
    ])
    def test_formats(self, numbers, expected):
        assert collapse_ranges(numbers) == expected


class TestStoredNumbering:
    """Mature residue numbers to chain positions, read off the file rather than computed."""

    def test_maps_through_the_file_offset(self, numbering):
        # mature 6 is the FIRST stored residue, so position 1 — not position 6
        assert numbering.ranges("A", [6], expect="A") == "1"
        assert numbering.ranges("A", [11], expect="G") == "6"

    def test_contiguous_mature_numbers_stay_contiguous(self, numbering):
        assert numbering.ranges("A", [6, 7, 8], expect="ACD") == "1..3"

    def test_refuses_when_the_residues_are_not_what_was_expected(self, numbering):
        """This is the guard: a shifted spec spells the wrong residues, and must not be emitted."""
        with pytest.raises(ValueError, match="numbering guard failed"):
            numbering.ranges("A", [6, 7, 8], expect="DEF")

    def test_the_error_names_both_sequences(self, numbering):
        with pytest.raises(ValueError) as caught:
            numbering.ranges("A", [6], expect="Z")
        assert "'A'" in str(caught.value) and "'Z'" in str(caught.value)

    def test_anti_hotspots_may_skip_the_check(self, numbering):
        """Omitting `expect` is allowed only where a shift costs a worse spec, not a wrong one."""
        assert numbering.ranges("A", [6, 7]) == "1..2"

    def test_sequence_reads_residues_by_mature_number(self, numbering):
        assert numbering.sequence("A", [6, 8, 10]) == "ADF"


class TestCofoldSpec:
    """The Boltz spec text."""

    def test_binder_is_chain_x_and_target_fills_a_to_c(self):
        spec = cofold_spec("TARGET", binder="BINDER")
        assert spec.startswith("version: 1\nsequences:\n")
        assert spec.count("msa: empty") == 4
        for chain in ("X", "A", "B", "C"):
            assert f"id: {chain}" in spec

    def test_target_only_spec_has_no_binder(self):
        spec = cofold_spec("TARGET")
        assert "id: X" not in spec and spec.count("msa: empty") == 3

    def test_protomer_count_is_honoured(self):
        assert cofold_spec("TARGET", protomers=1).count("msa: empty") == 1

    def test_chain_block_is_valid_yaml_shape(self):
        assert chain_block("A", "SEQ") == (
            "  - protein:\n      id: A\n      sequence: SEQ\n      msa: empty")


class TestHelixBundle:
    """The secondary-structure string BoltzGen is given for a helix-biased screen."""

    @pytest.mark.parametrize("n", [60, 70, 95, 110])
    def test_pattern_is_exactly_the_requested_length(self, n):
        """BoltzGen rejects a pattern longer than the sampled sequence, so this cannot drift."""
        assert len(helix_bundle(n)) == n

    def test_is_three_helices_separated_by_two_loops(self):
        """Uniform 'h' is not a designable topology -- the loops are the point."""
        runs = ["".join(g) for _, g in itertools.groupby(helix_bundle(60))]
        assert [r[0] for r in runs] == ["h", "l", "h", "l", "h"]
        assert [len(r) for r in runs if r[0] == "l"] == [4, 4]

    def test_sized_to_the_minimum_of_a_range_not_the_maximum(self):
        """A pattern cut to 95 would be rejected on every sample the range draws shorter."""
        assert len(secondary_structure("60..95")) == 60
