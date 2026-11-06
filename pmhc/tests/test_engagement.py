"""The peptide-vs-groove split, the chain-role partition, and the contact test underneath them."""

from __future__ import annotations

import gemmi
import numpy as np
import pytest
from common import PEPTIDE_G12D, PEPTIDE_WT, chain_roles, variant_of
from engagement import CONTACT, Engagement, contact_pairs

EPITOPES = (PEPTIDE_G12D, PEPTIDE_WT)

THREE_LETTER = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS",
    "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN",
    "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def structure_of(chains: dict[str, str]) -> gemmi.Structure:
    """A minimal multi-chain pose: one CA per residue, chains spaced far enough not to touch."""
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    for offset, (name, sequence) in enumerate(chains.items()):
        chain = gemmi.Chain(name)
        for index, letter in enumerate(sequence, start=1):
            residue = gemmi.Residue()
            residue.name = THREE_LETTER[letter]
            residue.seqid = gemmi.SeqId(index, " ")
            atom = gemmi.Atom()
            atom.name, atom.element = "CA", gemmi.Element("C")
            atom.pos = gemmi.Position(100.0 * offset, 3.8 * index, 0.0)
            residue.add_atom(atom)
            chain.add_residue(residue)
        structure[0].add_chain(chain)
    return structure


def row(peptide: int, hla: int) -> Engagement:
    return Engagement(
        design="d",
        variant="bulge_core",
        binder_length=80,
        peptide_contacts=peptide,
        hla_contacts=hla,
        hotspot_contacts=peptide,
        mutation_contacts=0,
        peptide_residues=3,
        hotspots_touched=3,
        hotspots_nominated=3,
    )


class TestPeptideShare:
    @pytest.mark.parametrize(
        "peptide,hla,expected",
        [
            # 7OW6's receptor, measured by this code: the one binder known to discriminate
            # these peptides puts a fifth of its interface on the peptide, so any gate
            # demanding most of it would reject the only molecule known to work.
            (55, 220, 0.2),
            # Same share at twice the interface: a raw contact count could not tell a
            # peptide-reader from a groove-gripper that merely binds harder.
            (110, 440, 0.2),
            (0, 300, 0.0),
            (50, 50, 0.5),
            (0, 0, 0.0),
        ],
    )
    def test_share(self, peptide, hla, expected):
        assert row(peptide, hla).peptide_share == pytest.approx(expected)

    def test_the_cutoff_travels_with_the_measurement(self):
        """A share measured at another cutoff is not comparable to the calibrated 0.200 bar."""
        assert row(55, 220).row["contact_cutoff"] == CONTACT


class TestContactPairs:
    @pytest.mark.parametrize(
        "distance,touching",
        [(CONTACT - 0.1, True), (CONTACT, True), (CONTACT + 0.1, False)],
    )
    def test_the_cutoff_is_inclusive(self, distance, touching):
        binder = np.array([[0.0, 0.0, 0.0]])
        target = np.array([[0.0, 0.0, distance]])
        assert bool(contact_pairs(binder, target)[0, 0]) is touching

    @pytest.mark.parametrize("empty_first", [True, False])
    def test_an_empty_side_gives_an_empty_matrix_rather_than_raising(self, empty_first):
        atoms = np.array([[0.0, 0.0, 0.0]])
        blank = np.empty((0, 3))
        pairs = contact_pairs(blank, atoms) if empty_first else contact_pairs(atoms, blank)
        assert pairs.sum() == 0
        assert pairs.shape == ((0, 1) if empty_first else (1, 0))

    def test_matches_the_explicit_difference_on_realistic_shapes(self):
        """The Gram identity is an optimisation, so it has to agree with the form it replaced
        pair for pair -- including at coordinates offset from the origin, where it loses the
        most precision."""
        rng = np.random.default_rng(0)
        for offset in (0.0, 100.0):
            binder = rng.normal(offset, 12.0, size=(220, 3))
            target = rng.normal(offset, 12.0, size=(900, 3))
            gaps = binder[:, None, :] - target[None, :, :]
            explicit = (gaps**2).sum(axis=2) <= CONTACT**2
            assert np.array_equal(contact_pairs(binder, target), explicit)


class TestChainRoles:
    def test_partitions_a_design_pose_by_sequence(self):
        structure = structure_of({"A": "GGGWWWGGG", "B": "KKKLLLKKK", "C": PEPTIDE_G12D})
        roles = chain_roles(
            structure, EPITOPES, targets=frozenset({"KKKLLLKKK", PEPTIDE_G12D})
        )
        assert (roles.design, roles.peptide, roles.receptor) == (("A",), "C", ("B",))

    def test_a_multi_chain_binder_is_named_not_inferred(self):
        """The 7OW6 control's binder is two chains, and it cannot be resolved the way a design
        is: the receptor sits in the novelty references on purpose, so nothing is left over."""
        structure = structure_of(
            {"D": "WWWGGG", "E": "YYYGGG", "A": "KKKLLL", "C": PEPTIDE_G12D}
        )
        roles = chain_roles(structure, EPITOPES, design=("D", "E"))
        assert (roles.design, roles.peptide, roles.receptor) == (("D", "E"), "C", ("A",))

    def test_the_wild_type_epitope_is_found_too(self):
        structure = structure_of({"A": "GGGWWW", "C": PEPTIDE_WT})
        assert chain_roles(structure, EPITOPES, design=("A",)).peptide == "C"

    @pytest.mark.parametrize(
        "chains,design,message",
        [
            ({"A": "GGGWWW", "B": "KKKLLL"}, ("A",), "epitope"),
            ({"A": "GGGWWW", "C": PEPTIDE_G12D}, ("Z",), "no such chain"),
            ({"C": PEPTIDE_G12D, "A": "GGGWWW"}, ("C",), "both the epitope"),
        ],
    )
    def test_refuses_a_partition_it_cannot_make(self, chains, design, message):
        with pytest.raises(ValueError, match=message):
            chain_roles(structure_of(chains), EPITOPES, design=design)

    def test_needs_either_named_chains_or_target_sequences(self):
        with pytest.raises(ValueError, match="either"):
            chain_roles(structure_of({"C": PEPTIDE_G12D}), EPITOPES)


class TestVariantOf:
    @pytest.mark.parametrize(
        "design,expected",
        [
            ("bulge_core_13", "bulge_core"),
            ("bulge_wide_07", "bulge_wide"),
            ("mut_proximal_59", "mut_proximal"),
            ("hub_39", "hub"),
        ],
    )
    def test_resolves_the_variant_from_the_design_name(self, design, expected):
        assert variant_of(design) == expected

    def test_rejects_a_name_whose_prefix_is_not_a_variant(self):
        """A rank-prefixed name is the shape that breaks this, so it must fail loudly."""
        with pytest.raises(ValueError):
            variant_of("rank01_bulge_core_25")
