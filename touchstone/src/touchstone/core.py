"""The contract between the two halves of touchstone.

A `Generator` (any of them) produces `BinderDesign`s; a `Verifier` (the asset)
consumes them and returns a `Verdict`. The verifier never sees *how* a design was
made — only the `CoordinationSite` — which is what makes it generator-blind.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import combinations
from typing import Protocol, runtime_checkable

import numpy as np


def element_symbol(metal: str) -> str:
    """ASE element symbol for a site-metal label: 'Ni2+' -> 'Ni'."""
    return "".join(c for c in metal if c.isalpha())


def oxidation_state(metal: str) -> int:
    """Formal charge from a site-metal label: 'Ni2+' -> 2, 'Fe3+' -> 3, 'Cu+' -> 1."""
    m = re.search(r"(\d*)\s*([+-])", metal)  # bare sign (e.g. 'Cu+') ⇒ magnitude 1
    if not m:
        raise ValueError(f"no oxidation state in {metal!r}")
    return int(m.group(1) or 1) * (1 if m.group(2) == "+" else -1)


# Programming errors — a typo, wrong attribute, wrong call, or unimplemented stub. A
# verifier must let these surface, not launder them into a benign "defer" that reads
# like a missing backend.
_PROGRAMMING_ERRORS = (AttributeError, TypeError, NameError, NotImplementedError)


def reraise_if_bug(e: BaseException) -> None:
    """Re-raise programming errors so they surface; call inside `except Exception as e`
    before deferring, so only genuine runtime/backend failures degrade to a defer."""
    if isinstance(e, _PROGRAMMING_ERRORS):
        raise e


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

    @property
    def is_empty(self) -> bool:
        return self.coordination_number == 0

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
    source: str | None = None  # path the design was loaded from, if any

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

    @property
    def label(self) -> str:
        """Single-word verdict for display/grouping: defer / trust / weak."""
        return "defer" if self.ood else ("trust" if self.trust else "weak")

    @classmethod
    def defer(cls, reason: str, score: float = 0.0) -> "Verdict":
        """An off-manifold verdict — not trusted, flagged for review. Owns the
        '— defer' reason suffix so every verifier signals it the same way."""
        return cls(score, trust=False, ood=True, reason=f"{reason} — defer")


@runtime_checkable
class Generator(Protocol):
    def design(self, target: str, n: int = ...) -> list[BinderDesign]: ...


@runtime_checkable
class Verifier(Protocol):
    def verify(self, design: BinderDesign) -> Verdict: ...
