"""Chain classification, mature numbering and buried area, on a synthetic complex.

A synthetic rather than a real complex so the geometry is known: the binder is placed against
one end of the target, so burial MUST be local, and the target chain is an exact slice of human
TNF-alpha at a known offset, so the mature numbering it should recover is known too.

These three are the inputs to every interface verdict, and `verify` has no other coverage.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.PDB import PDBIO, StructureBuilder

from tnf_binder.structure import Complex, chain_sequence, contains_tnf
from tnf_binder.target import HUMAN

TARGET_START = 100  # 0-based slice into HUMAN, so the first target residue is mature 101
TARGET_LENGTH = 30
BINDER = "LEEKLRE"
RESIDUE_RISE = 3.5  # A between consecutive residues along y
BACKBONE = (("N", 0.0), ("CA", 1.2), ("C", 2.4), ("O", 3.0))


def build_complex(path: Path, target_seq: str, binder_seq: str, gap: float) -> Path:
    """Two extended chains side by side, `gap` Angstrom apart along x."""
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure("synthetic")
    builder.init_model(0)
    for chain_id, sequence, x in (("A", target_seq, 0.0), ("B", binder_seq, gap)):
        builder.init_chain(chain_id)
        builder.init_seg(" ")
        for number, residue in enumerate(sequence, start=1):
            builder.init_residue(protein_letters_1to3[residue].upper(), " ", number, " ")
            for name, z in BACKBONE:
                builder.init_atom(name, (x, number * RESIDUE_RISE, z), 0.0, 1.0, " ",
                                  name, element=name[0])
    io = PDBIO()
    io.set_structure(builder.get_structure())
    io.save(str(path))
    return path


@pytest.fixture
def complex_(tmp_path) -> Complex:
    target = HUMAN[TARGET_START:TARGET_START + TARGET_LENGTH]
    return Complex.from_file(build_complex(tmp_path / "c.pdb", target, BINDER, gap=4.0))


class TestChainClassification:
    """Which chain is the target is decided by TNF substring containment, not identity."""

    def test_target_and_binder_are_separated(self, complex_):
        assert complex_.is_resolved
        assert [c.id for c in complex_.target_chains] == ["A"]
        assert complex_.binder_sequence == BINDER

    def test_a_short_binder_is_never_mistaken_for_the_target(self):
        assert contains_tnf(HUMAN[TARGET_START:TARGET_START + TARGET_LENGTH])
        assert not contains_tnf(BINDER)

    def test_chain_sequence_round_trips(self, complex_):
        target = HUMAN[TARGET_START:TARGET_START + TARGET_LENGTH]
        assert chain_sequence(complex_.target_chains[0]) == target


class TestMatureNumbering:
    """Stored numbering is not trusted; positions come from aligning to human TNF-alpha."""

    def test_recovers_the_offset_the_slice_was_taken_at(self, complex_):
        numbering = complex_.mature_numbering()
        positions = sorted(numbering.values())
        assert positions == list(range(TARGET_START + 1,
                                      TARGET_START + 1 + TARGET_LENGTH))

    def test_numbering_is_not_the_stored_residue_number(self, complex_):
        """The file numbers these residues 1..30; mature numbering must not echo that."""
        stored = sorted(r.id[1] for r in complex_.target_chains[0])
        assert stored[0] == 1
        assert min(complex_.mature_numbering().values()) == TARGET_START + 1


class TestBuriedByPosition:
    """Buried area, keyed by mature position."""

    def test_burial_is_local_to_where_the_binder_sits(self, complex_):
        areas = complex_.buried_by_position()
        assert areas, "a binder 4 A away should bury something"
        # the 7-residue binder spans the first ~7 target residues, i.e. mature 101-108
        assert max(areas) <= TARGET_START + len(BINDER) + 2
        assert all(area > 0 for area in areas.values())

    def test_the_far_end_of_the_target_is_untouched(self, complex_):
        areas = complex_.buried_by_position()
        far_end = TARGET_START + TARGET_LENGTH
        assert far_end not in areas

    def test_a_distant_binder_buries_nothing(self, tmp_path):
        target = HUMAN[TARGET_START:TARGET_START + TARGET_LENGTH]
        far = Complex.from_file(build_complex(tmp_path / "far.pdb", target, BINDER, gap=40.0))
        assert far.buried_by_position() == {}
        assert far.clashes() == 0

    def test_an_overlapping_binder_clashes(self, tmp_path):
        target = HUMAN[TARGET_START:TARGET_START + TARGET_LENGTH]
        overlapped = Complex.from_file(
            build_complex(tmp_path / "clash.pdb", target, BINDER, gap=0.5))
        assert overlapped.clashes() > 0
