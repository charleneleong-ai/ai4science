"""Molecular formula → monoisotopic mass.

Pool masses must come from the formula, not from an observed precursor m/z. Deriving them
from the observation makes the true candidate's mass float-identical to the query's, so a
mass-error ranker recovers it by measurement identity rather than chemistry — measured at
100% of truths sitting at exactly 0.0000 ppm, against a 1.69 ppm decoy median.

From the formula, every isomer of a formula shares one mass, so mass error carries no
signal between them. That is the honest situation the harness has to reproduce.
"""

from __future__ import annotations

import re

# most-abundant-isotope masses
MONOISOTOPIC: dict[str, float] = {
    "H": 1.00782503, "D": 2.01410178, "C": 12.0, "N": 14.00307400, "O": 15.99491462,
    "F": 18.99840320, "Na": 22.98976928, "Si": 27.97692653, "P": 30.97376151,
    "S": 31.97207069, "Cl": 34.96885271, "K": 38.96370649, "Ca": 39.96259086,
    "Fe": 55.93494200, "Co": 58.93319400, "Ni": 57.93534241, "Cu": 62.92959772,
    "Zn": 63.92914201, "As": 74.92159460, "Se": 79.91651960, "Br": 78.91833760,
    "I": 126.90446800, "B": 11.00930550, "Mg": 23.98504170, "Mn": 54.93804391,
    "Hg": 201.97064340, "Pt": 194.96479170, "Sn": 119.90220160,
}

TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> dict[str, int] | None:
    """Element counts, or None if the formula contains anything unrecognised.

    None rather than a partial parse: silently dropping an element understates the mass and
    moves the candidate into a different retrieval window.
    """
    if not formula:
        return None
    cleaned = formula.strip().replace(" ", "")
    # Strip trailing charge signs only. A digit before the sign is ambiguous — in "C6H12O6+"
    # the 6 is an element count, not a charge magnitude — so leave digits alone rather than
    # risk silently halving a mass.
    cleaned = cleaned.rstrip("+-")
    if not cleaned:
        return None
    counts: dict[str, int] = {}
    pos = 0
    for match in TOKEN.finditer(cleaned):
        if match.start() != pos:  # a gap means an unparsed character
            return None
        pos = match.end()
        element, digits = match.group(1), match.group(2)
        if element not in MONOISOTOPIC:
            return None
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    if pos != len(cleaned):
        return None
    return counts or None


def formula_mass(formula: str) -> float | None:
    """Neutral monoisotopic mass of a molecular formula, or None if unparseable."""
    counts = parse_formula(formula)
    if counts is None:
        return None
    return sum(MONOISOTOPIC[el] * n for el, n in counts.items())
