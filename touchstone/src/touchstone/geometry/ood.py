"""Extreme-condition perturbation — the OOD probe.

Real critical-mineral recovery happens in hot, acidic, saline leachate. Protonation
weakens metal–ligand bonds, lengthening them. `under_leachate` applies that proxy so
the verifier can answer: does this binder stay on the reference manifold under the
conditions that break it, or does it go out-of-distribution (→ defer)?
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from ..core import BinderDesign, CoordinationSite


def acidic_leachate(site: CoordinationSite, bond_stretch: float = 0.5) -> CoordinationSite:
    """Stretch every metal–ligand bond radially outward by `bond_stretch` Å."""
    direction = site.ligand_xyz - site.metal_xyz
    unit = direction / np.linalg.norm(direction, axis=1, keepdims=True)
    return CoordinationSite(
        metal=site.metal,
        metal_xyz=site.metal_xyz,
        ligand_xyz=site.ligand_xyz + unit * bond_stretch,
        ligand_elems=site.ligand_elems,
    )


def under_leachate(design: BinderDesign, bond_stretch: float = 0.5) -> BinderDesign:
    """A copy of `design` with its site perturbed toward extreme leachate."""
    return replace(design, site=acidic_leachate(design.site, bond_stretch))
