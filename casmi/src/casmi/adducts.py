"""Precursor m/z → neutral monoisotopic mass, per adduct.

Retrieval windows are computed in neutral-mass space, so a wrong shift here moves every
candidate set. Values are monoisotopic and electron-corrected where the charge demands it.
"""

from __future__ import annotations

PROTON = 1.007276

# adduct → (multiplier on M, m/z offset). neutral M = (mz - offset) / multiplier
ADDUCTS: dict[str, tuple[int, float]] = {
    "[M+H]+": (1, PROTON),
    "[M-H]-": (1, -PROTON),
    "[M+Na]+": (1, 22.989218),
    "[M+NH4]+": (1, 18.033823),
    "[M+K]+": (1, 38.963158),
    # formic-acid adduct: M + HCOOH - H
    "[M+CH2O2-H]-": (1, 46.005480 - PROTON),
    "[M+Cl]-": (1, 34.969402),
    # dimers: the precursor carries two copies of M
    "[2M+H]+": (2, PROTON),
    "[2M+Na]+": (2, 22.989218),
}


def neutral_mass(precursor_mz: float, adduct: str) -> float | None:
    """The neutral monoisotopic mass implied by a precursor, or None for an unmapped adduct.

    None rather than a guess: an unknown adduct silently mapped to [M+H]+ would shift the
    retrieval window and quietly drop the true candidate.
    """
    entry = ADDUCTS.get(adduct)
    if entry is None:
        return None
    multiplier, offset = entry
    return (precursor_mz - offset) / multiplier


def ppm_window(mass: float, ppm: float) -> tuple[float, float]:
    tol = mass * ppm * 1e-6
    return mass - tol, mass + tol
