"""touchstone — a generator-agnostic geometry verifier for designed metal binders.

A touchstone was historically used to assay precious metals: the independent
reference standard against which a sample is judged. Same job here.
"""

from .core import BinderDesign, CoordinationSite, Generator, Verdict, Verifier
from .generators import BoltzGenAdapter, MockGenerator, RFdiffusionAdapter, octahedral_site
from .geometry.ood import under_leachate
from .geometry.parse import coordination_site_from_pdb
from .geometry.reference import MetalGeometry, MockReference, PDBReference, ReferenceDistribution
from .geometry.verifier import GeometryVerifier
from .pipeline import design_and_rank, rank

__all__ = [
    "BinderDesign",
    "CoordinationSite",
    "Generator",
    "Verdict",
    "Verifier",
    "MockGenerator",
    "RFdiffusionAdapter",
    "BoltzGenAdapter",
    "octahedral_site",
    "under_leachate",
    "coordination_site_from_pdb",
    "MetalGeometry",
    "MockReference",
    "PDBReference",
    "ReferenceDistribution",
    "GeometryVerifier",
    "design_and_rank",
    "rank",
]
