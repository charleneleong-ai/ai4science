"""Geometry validation against the CSD via Mogul (CCDC) — the gold-standard
counterpart to GeometryVerifier's hand-rolled z-score.

Mogul checks each metal–donor bond against the Cambridge Structural Database and
reports how unusual the length is (z-score) and how many CSD hits backed the
comparison. This verifier trusts a site whose bonds are all within `trust_z`, defers
a fragment with too few CSD hits to judge (or when Mogul can't run), and flags an
off-distribution bond. It uses CCDC's *validation engine*, not just CSD data — so it
is license-gated (CSD Python API); the analyser is pluggable (inject one for tests,
the default lazily drives `ccdc.conformer.GeometryAnalyser`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..core import BinderDesign, CoordinationSite, Verdict, element_symbol


@dataclass
class MogulFragment:
    """One Mogul-checked geometry (a metal–donor bond) against the CSD."""

    label: str  # e.g. "Ni-N"
    value: float  # observed length, Å
    z_score: float  # (value − CSD mean) / CSD std
    nhits: int  # CSD sample size behind the comparison


def _ccdc_analyse(site: CoordinationSite, metal: str) -> list[MogulFragment]:
    """Drive CCDC's Mogul (GeometryAnalyser) over the metal–donor bonds. Requires the
    CSD Python API (`ccdc`) + a licence — the exact result attributes should be
    confirmed against the installed API version before quantitative use.
    """
    from ccdc.conformer import GeometryAnalyser
    from ccdc.molecule import Atom, Bond, Molecule

    mol = Molecule("site")
    m_atom = mol.add_atom(Atom(metal, coordinates=tuple(map(float, site.metal_xyz))))
    for elem, xyz in zip(site.ligand_elems, site.ligand_xyz):
        a = mol.add_atom(Atom(elem, coordinates=tuple(map(float, xyz))))
        mol.add_bond(Bond.BondType(1), m_atom, a)

    analysed = GeometryAnalyser().analyse_molecule(mol)
    return [
        MogulFragment(label="-".join(b.atom_labels), value=b.value, z_score=b.z_score, nhits=b.nhits)
        for b in analysed.analysed_bonds
    ]


class MogulVerifier:
    """Trusts a site whose metal–donor bonds are all geometrically normal vs the CSD
    (Mogul |z| ≤ `trust_z`); defers when CSD support is too thin to judge or Mogul
    can't run; flags a bond beyond `ood_z` as off-distribution.
    """

    def __init__(
        self,
        analyse: Callable[[CoordinationSite, str], list[MogulFragment]] | None = None,
        *,
        trust_z: float = 2.0,
        ood_z: float = 5.0,
        min_hits: int = 15,
    ):
        self._analyse = analyse or _ccdc_analyse
        self.trust_z = trust_z
        self.ood_z = ood_z
        self.min_hits = min_hits

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            frags = self._analyse(design.site, element_symbol(design.site.metal))
        except Exception as e:  # no licence / Mogul failure ⇒ can't validate
            return Verdict.defer(f"Mogul analysis unavailable: {type(e).__name__}")
        if not frags:
            return Verdict.defer("no Mogul-matched fragments")
        if (fewest := min(f.nhits for f in frags)) < self.min_hits:
            return Verdict.defer(f"insufficient CSD support ({fewest} < {self.min_hits} hits)")

        max_z = max(abs(f.z_score) for f in frags)
        score = float(np.exp(-0.5 * (max_z / self.trust_z) ** 2))
        reason = f"Mogul max |z| {max_z:.1f}σ over {len(frags)} bonds ({fewest}+ CSD hits)"
        if max_z > self.ood_z:
            return Verdict.defer(reason, score=score)
        return Verdict(score, trust=max_z <= self.trust_z, ood=False, reason=reason)
