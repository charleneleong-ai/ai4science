"""The contract between the two halves of touchstone.

A `Generator` (any of them) produces `BinderDesign`s; a `Verifier` (the asset)
consumes them and returns a `Verdict`. The verifier never sees *how* a design was
made — only the `CoordinationSite` — which is what makes it generator-blind.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass
class CoordinationSite:
    """A metal centre and the atoms coordinating it, in Angstrom."""

    metal: str  # e.g. "Ni2+"
    metal_xyz: np.ndarray  # (3,)
    ligand_xyz: np.ndarray  # (n, 3)
    ligand_elems: tuple[str, ...]  # coordinating-atom elements, e.g. ("N", "O", ...)

    @property
    def coordination_number(self) -> int:
        return len(self.ligand_xyz)

    def bond_lengths(self) -> np.ndarray:
        """Metal–ligand distances."""
        return np.linalg.norm(self.ligand_xyz - self.metal_xyz, axis=1)

    def bond_angles(self) -> np.ndarray:
        """All ligand–metal–ligand angles, in degrees."""
        v = self.ligand_xyz - self.metal_xyz
        v = v / np.linalg.norm(v, axis=1, keepdims=True)
        cos = [np.clip(np.dot(v[i], v[j]), -1.0, 1.0) for i, j in combinations(range(len(v)), 2)]
        return np.degrees(np.arccos(cos)) if cos else np.empty(0)


@dataclass
class BinderDesign:
    """What a generator emits and the verifier consumes — the only coupling point."""

    sequence: str
    site: CoordinationSite
    generator: str  # which generator produced this (provenance, not trusted)
    generator_confidence: float  # the generator's *own* score — recorded, never trusted

    @property
    def target_metal(self) -> str:
        return self.site.metal


@dataclass
class Verdict:
    """The verifier's call on one design. Higher `score` = more geometrically plausible."""

    score: float
    trust: bool  # plausible geometry, in-distribution
    ood: bool  # site sits off the reference manifold (e.g. extreme-leachate input)
    reason: str


@runtime_checkable
class Generator(Protocol):
    def design(self, target: str, n: int = ...) -> list[BinderDesign]: ...


@runtime_checkable
class Verifier(Protocol):
    def verify(self, design: BinderDesign) -> Verdict: ...
