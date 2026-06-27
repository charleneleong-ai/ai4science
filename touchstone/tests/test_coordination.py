import numpy as np

from touchstone import BinderDesign, CoordinationGeometryVerifier, CoordinationSymmetryVerifier, octahedral_site
from touchstone.core import CoordinationSite


def _design(site: CoordinationSite) -> BinderDesign:
    return BinderDesign("", site, generator="t", generator_confidence=0.5)


def _site(xyz: list[list[float]], elems: tuple[str, ...]) -> CoordinationSite:
    return CoordinationSite("Ni2+", np.zeros(3), np.array(xyz), elems)


class TestCoordinationSymmetry:
    """nVECSUM: is the metal enclosed, or are the donors all on one side?"""

    def test_octahedron_is_balanced(self):
        v = CoordinationSymmetryVerifier().verify(_design(octahedral_site("Ni2+")))
        assert v.trust and v.metrics["nvecsum"] < 0.05  # ±x,±y,±z unit vectors cancel

    def test_one_sided_coordination_defers(self):
        # three donors crammed onto the +x hemisphere ⇒ metal half-exposed
        v = CoordinationSymmetryVerifier().verify(_design(_site([[2, 0, 0], [1.8, 0.8, 0], [1.8, -0.8, 0]], ("N", "N", "O"))))
        assert v.ood and not v.trust and v.metrics["nvecsum"] > 0.6


class TestCoordinationGeometry:
    """Polyhedron shape: right distances + count, but the right arrangement?"""

    def test_octahedron_matches_ideal(self):
        v = CoordinationGeometryVerifier().verify(_design(octahedral_site("Ni2+")))
        assert v.trust and v.metrics["angle_rmsd_deg"] < 1.0

    def test_tetrahedron_matches_ideal(self):
        # four tetrahedron vertices ⇒ ~109.5° angles ⇒ the tetrahedral isomer wins over square-planar
        t = 2.0 / np.sqrt(3)
        v = CoordinationGeometryVerifier().verify(
            _design(_site([[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]], ("N", "N", "O", "O")))
        )
        assert v.trust and v.metrics["angle_rmsd_deg"] < 2.0

    def test_collapsed_geometry_defers(self):
        # four donors squeezed into a narrow cone — nothing like a tetrahedron/square-plane
        v = CoordinationGeometryVerifier().verify(
            _design(_site([[2, 0, 0], [1.9, 0.4, 0], [1.9, 0, 0.4], [1.9, 0.4, 0.4]], ("N", "N", "O", "O")))
        )
        assert v.ood and v.metrics["angle_rmsd_deg"] > 40

    def test_no_ideal_reference_for_high_cn_defers(self):
        xyz = [[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2], [1.4, 1.4, 0]]  # CN7
        v = CoordinationGeometryVerifier().verify(_design(_site(xyz, ("N",) * 7)))
        assert v.ood and "CN=7" in v.reason
