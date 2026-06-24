from pathlib import Path

import numpy as np
import pytest

from touchstone import coordination_site_from_pdb

FIXTURE = Path(__file__).parent / "fixtures" / "rfaa_nickel_sample_0.pdb"


def _atom(serial, name, resname, x, y, z, element, rec="HETATM"):
    return (
        f"{rec:<6}{serial:>5} {name:<4} {resname:>3} A{1:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
    )


def _write_site(tmp_path, donors, far=(), nondonor=True):
    """NI at origin + octahedral `donors` at 2.1 Å, optional far donors + a carbon."""
    lines = [_atom(1, "NI", "NI", 0, 0, 0, "NI")]
    for i, (el, vec) in enumerate(donors, start=2):
        x, y, z = np.array(vec) * 2.1
        lines.append(_atom(i, el, "HIS", x, y, z, el, rec="ATOM"))
    for j, (el, d) in enumerate(far, start=100):
        lines.append(_atom(j, el, "ASP", d, 0, 0, el, rec="ATOM"))
    if nondonor:
        lines.append(_atom(200, "C", "ALA", 2.1, 0, 0, "C", rec="ATOM"))  # not a donor
    p = tmp_path / "design.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


OCTAHEDRON = [
    ("N", (1, 0, 0)), ("N", (-1, 0, 0)), ("O", (0, 1, 0)),
    ("O", (0, -1, 0)), ("N", (0, 0, 1)), ("O", (0, 0, -1)),
]


class TestCoordinationSiteFromPDB:
    def test_extracts_metal_and_shell(self, tmp_path):
        site = coordination_site_from_pdb(_write_site(tmp_path, OCTAHEDRON), "NI", "Ni2+")
        assert site.metal == "Ni2+"
        assert site.coordination_number == 6
        assert np.allclose(site.bond_lengths(), 2.1)

    def test_excludes_far_donors_and_nondonors(self, tmp_path):
        # a far O at 5 Å and a C at 2.1 Å must not enter the first shell
        pdb = _write_site(tmp_path, OCTAHEDRON, far=[("O", 5.0)], nondonor=True)
        site = coordination_site_from_pdb(pdb, "NI", "Ni2+", cutoff=2.8)
        assert site.coordination_number == 6
        assert "C" not in site.ligand_elems

    def test_cutoff_controls_shell_size(self, tmp_path):
        pdb = _write_site(tmp_path, OCTAHEDRON, far=[("O", 2.5)])
        assert coordination_site_from_pdb(pdb, "NI", "Ni2+", cutoff=2.3).coordination_number == 6
        assert coordination_site_from_pdb(pdb, "NI", "Ni2+", cutoff=2.8).coordination_number == 7

    def test_excludes_placeholder_atoms_on_the_metal(self, tmp_path):
        # RFdiffusionAA parks unplaced sidechain atoms at the origin, coincident with
        # a recentred metal — these must not count as 0 Å bonds.
        lines = [_atom(1, "NI", "NI", 0, 0, 0, "NI")]
        lines.append(_atom(2, "OG", "SER", 0, 0, 0, "O", rec="ATOM"))  # unplaced placeholder
        for i, (el, vec) in enumerate(OCTAHEDRON, start=3):
            x, y, z = np.array(vec) * 2.1
            lines.append(_atom(i, el, "HIS", x, y, z, el, rec="ATOM"))
        p = tmp_path / "placeholder.pdb"
        p.write_text("\n".join(lines) + "\nEND\n")
        site = coordination_site_from_pdb(p, "NI", "Ni2+")
        assert site.coordination_number == 6
        assert (site.bond_lengths() > 0.5).all()

    def test_missing_metal_raises(self, tmp_path):
        p = tmp_path / "no_metal.pdb"
        p.write_text(_atom(1, "N", "HIS", 1, 1, 1, "N", rec="ATOM") + "\nEND\n")
        with pytest.raises(ValueError, match="no 'NI'"):
            coordination_site_from_pdb(p, "NI", "Ni2+")


class TestRealRFdiffusionOutput:
    """Parses an actual RFdiffusionAA nickel design (backbone-only, blank element
    column, sidechain placeholders parked on the metal)."""

    def test_real_nickel_design(self):
        site = coordination_site_from_pdb(FIXTURE, "NI", "Ni2+")
        assert site.coordination_number == 2  # two backbone N donors near the Ni
        assert set(site.ligand_elems) == {"N"}
        assert (site.bond_lengths() > 1.0).all()  # placeholders at the origin excluded
