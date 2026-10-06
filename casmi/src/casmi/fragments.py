"""Rank candidates by how much of the observed spectrum their structure can explain.

The harness measured the real problem: a 5 ppm window holds ~46 candidates of which ~29
share the truth's exact formula. Mass accuracy cannot separate those, so the ranker needs a
signal that depends on *connectivity*. Breaking bonds does.

Combinatorial fragmentation (the MetFrag idea) is deliberately the first attempt: no
training, no spectral library, no model weights, so it establishes whether connectivity
carries the signal before anything heavier is justified.

Ring bonds are breakable. A first version skipped them, and measured +0.006 over the
mass-only floor — because 68% of bonds in these natural products are ring bonds, so only
~19% of peak intensity could be explained at all. Opening a ring needs two cuts in it,
which is why `max_breaks` must be >= 2 to see ring fragments.

Known blind spot, pinned in tests: isomers differing only by branching cleave into pieces of
the *same* formula, so no amount of enumeration separates them.
"""

from __future__ import annotations

import bisect
from collections.abc import Sequence
from functools import lru_cache
from itertools import combinations

from rdkit import Chem, RDLogger

from .adducts import PROTON
from .formula import MONOISOTOPIC
from .harness import Query, Structure

RDLogger.DisableLog("rdApp.*")  # candidate SMILES come from a database; parse noise is expected

H_MASS = 1.00782503
# ions are observed, not neutrals, and fragments rearrange hydrogens — so a fragment may
# appear at +-1 H either side of its bare mass
H_SHIFTS = (-1, 0, 1)
MAX_BONDS = 60   # past this the pair enumeration stops paying for itself
MAX_BREAKS = 2


class MolGraph:
    """Atom masses plus a bond list — enough to take fragment masses by connectivity.

    Built once per structure instead of calling RDKit's fragmenter per bond combination,
    which is what made ring enumeration affordable: including ring bonds multiplies the
    number of combinations by roughly ten.
    """

    __slots__ = ("masses", "bonds", "total")

    def __init__(self, mol: Chem.Mol) -> None:
        self.masses: list[float] = [
            MONOISOTOPIC.get(a.GetSymbol(), 0.0) + a.GetTotalNumHs() * H_MASS
            for a in mol.GetAtoms()
        ]
        self.bonds: list[tuple[int, int]] = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
        ]
        self.total: float = sum(self.masses)


def component_masses(graph: MolGraph, cut: frozenset[int]) -> list[float]:
    """Masses of the pieces left when these bond indices are removed."""
    parent = list(range(len(graph.masses)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, (a, b) in enumerate(graph.bonds):
        if i in cut:
            continue
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    sums: dict[int, float] = {}
    for idx, mass in enumerate(graph.masses):
        root = find(idx)
        sums[root] = sums.get(root, 0.0) + mass
    # round before the caller sets them: summing in different orders yields values differing
    # in the last bits, which a float set keeps as distinct entries and bloats every lookup.
    # 1e-6 Da is four orders below the 0.01 Da match tolerance, so nothing is lost.
    return [round(v, 6) for v in sums.values()]


def fragment_mass_set(smiles: str, max_breaks: int = MAX_BREAKS,
                      include_rings: bool = True) -> frozenset[float]:
    """Neutral monoisotopic masses of every fragment reachable by cutting <= max_breaks
    bonds, plus the intact molecule.

    `include_rings=False` reproduces the first version, for comparison only — it cannot see
    ring fragments, which is most of the signal in a natural product.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return frozenset()
    graph = MolGraph(mol)
    masses = {round(graph.total, 6)}
    breakable = [
        i for i, _ in enumerate(graph.bonds)
        if include_rings or not mol.GetBondWithIdx(i).IsInRing()
    ]
    if len(breakable) > MAX_BONDS:
        return frozenset(masses)
    for n in range(1, max_breaks + 1):
        for combo in combinations(breakable, n):
            masses.update(component_masses(graph, frozenset(combo)))
    return frozenset(masses)


@lru_cache(maxsize=200_000)
def cached_fragments(smiles: str) -> frozenset[float]:
    """Candidates recur across queries, so enumerate each structure once."""
    return fragment_mass_set(smiles)


def explained_fraction(mzs: Sequence[float], intensities: Sequence[float],
                       frag_masses: frozenset[float], positive: bool, tol_da: float = 0.01) -> float:
    """Share of total peak intensity explainable as a fragment of this structure.

    Intensity-weighted rather than a peak count: explaining one base peak says more about a
    structure than explaining several near-noise peaks.
    """
    if not mzs or not frag_masses:
        return 0.0
    charge = PROTON if positive else -PROTON
    total = sum(intensities) or 1.0
    hit = 0.0
    sorted_frags = sorted(frag_masses)
    for mz, inten in zip(mzs, intensities, strict=True):
        if any_within(sorted_frags, mz - charge, tol_da):
            hit += inten
    return hit / total


def any_within(sorted_masses: Sequence[float], target: float, tol_da: float) -> bool:
    """Is some fragment mass within tolerance of the target, allowing +-1 H rearrangement?"""
    for shift in H_SHIFTS:
        want = target + shift * H_MASS
        i = bisect.bisect_left(sorted_masses, want - tol_da)
        if i < len(sorted_masses) and sorted_masses[i] <= want + tol_da:
            return True
    return False


def score_candidate(query: Query, structure: Structure) -> float:
    """Best explained fraction across the molecule's spectra.

    Max, not mean: a molecule's spectra span collision energies, and a low-energy spectrum
    that barely fragments should not penalise a structure the high-energy one supports.
    """
    frags = cached_fragments(structure.smiles)
    if not frags:
        return 0.0
    positive = query.adduct.endswith("+")
    return max(
        (explained_fraction(mzs, intens, frags, positive) for mzs, intens in query.spectra),
        default=0.0,
    )


def rank_by_fragments(query: Query, candidates: Sequence[Structure]) -> list[Structure]:
    """Most-explained first; mass error breaks ties so the ordering is never arbitrary."""
    mass = query.mass or 0.0
    return sorted(candidates, key=lambda s: (-score_candidate(query, s), abs(s.mass - mass)))
