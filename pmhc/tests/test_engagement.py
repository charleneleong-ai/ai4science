"""The peptide-vs-groove split, and the contact test underneath it."""

from __future__ import annotations

import numpy as np
import pytest
from engagement import CONTACT, Engagement, contact_pairs, variant_of


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
    def test_interface_size_does_not_change_the_share(self):
        """The point of a ratio: a bigger interface raises both counts, so a raw peptide-contact
        count cannot tell a peptide-reader from a groove-gripper that simply binds harder."""
        assert row(55, 220).peptide_share == row(110, 440).peptide_share

    @pytest.mark.parametrize(
        "peptide,hla,expected",
        [(55, 220, 0.2), (0, 300, 0.0), (50, 50, 0.5), (0, 0, 0.0)],
    )
    def test_share(self, peptide, hla, expected):
        assert row(peptide, hla).peptide_share == pytest.approx(expected)

    def test_the_known_binder_would_fail_an_intuitive_bar(self):
        """7OW6's receptor, measured by this code, is 55 peptide against 220 HLA. Any gate
        demanding most of the interface on the peptide rejects the one binder known to work."""
        assert row(55, 220).peptide_share < 0.5


class TestContactPairs:
    def test_counts_pairs_inside_the_cutoff_and_excludes_the_rest(self):
        binder = np.array([[0.0, 0.0, 0.0]])
        target = np.array([[0.0, 0.0, CONTACT - 0.1], [0.0, 0.0, CONTACT + 0.1]])
        assert contact_pairs(binder, target).tolist() == [[True, False]]

    def test_the_cutoff_is_inclusive(self):
        binder = np.array([[0.0, 0.0, 0.0]])
        assert contact_pairs(binder, np.array([[0.0, 0.0, CONTACT]])).all()

    @pytest.mark.parametrize("empty_first", [True, False])
    def test_an_empty_side_gives_an_empty_matrix_rather_than_raising(self, empty_first):
        atoms = np.array([[0.0, 0.0, 0.0]])
        blank = np.empty((0, 3))
        pairs = contact_pairs(blank, atoms) if empty_first else contact_pairs(atoms, blank)
        assert pairs.sum() == 0
        assert pairs.shape == ((0, 1) if empty_first else (1, 0))


class TestVariantOf:
    @pytest.mark.parametrize(
        "design,expected",
        [
            ("bulge_core_13", "bulge_core"),
            ("bulge_wide_07", "bulge_wide"),
            ("mut_proximal_59", "mut_proximal"),
            ("hub_39", "hub"),
            ("rank01_bulge_core_25", None),
        ],
    )
    def test_resolves_the_variant_from_the_design_name(self, design, expected):
        if expected is None:
            with pytest.raises(Exception):
                variant_of(design)
        else:
            assert variant_of(design) == expected
