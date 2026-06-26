"""MLIP physics tier: relax a metal-coordination cluster under a machine-learned
interatomic potential and judge whether the site holds.

Sharper than the semi-empirical xtb step — DFT-accuracy energies/forces at MD
speed — and backbone-pluggable: MACE-MP by default, UMA when its gated weights are
available, ASE EMT as a dependency-light fallback. The heavy backends are optional
and mutually exclusive (mace-torch and fairchem pin conflicting e3nn): install
`touchstone[mace]` *or* `touchstone[uma]`, in separate environments.

Operates on a protonated structure (same prep as the xtb flow — see
scripts/extract_cluster.py): extracts the metal-centred cluster, relaxes it, and
reports site stability plus a ranking-only metal interaction energy. ASE imports
are deferred so the module loads without the optional backends installed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core import BinderDesign, Verdict, element_symbol
from ..geometry.parse import DONOR_ELEMENTS


def make_backbone(name: str, device: str = "cuda"):
    """ASE calculator for a backbone. Heavy imports are deferred to call time."""
    if name == "mace_mp":
        from mace.calculators import mace_mp

        return mace_mp(model="medium", device=device, default_dtype="float64")
    if name == "uma":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        unit = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)
        return FAIRChemCalculator(unit, task_name="omol")
    if name == "emt":  # ASE-core, light — real relaxation for Ni/Cu/… smoke runs
        from ase.calculators.emt import EMT

        return EMT()
    raise ValueError(f"unknown MLIP backbone {name!r}")


@dataclass
class SiteRelaxation:
    """Metal-site geometry across an MLIP relaxation."""

    cn_before: int
    cn_after: int
    site_drift: float  # max displacement of metal + original first-shell atoms (Å)
    mean_bond: float  # mean metal–donor distance after relaxation (Å)
    interaction_energy: float | None  # E(complex) − E(apo) − E(metal), ranking-only (eV)

    @property
    def donors_lost(self) -> int:
        return max(0, self.cn_before - self.cn_after)


def _metal_index(symbols: list[str], metal: str) -> int:
    if metal not in symbols:
        raise ValueError(f"no {metal} atom in cluster")
    return symbols.index(metal)


def _shell(positions: np.ndarray, symbols: list[str], mi: int, cutoff: float) -> dict[int, float]:
    """{donor index: distance} for N/O/S within cutoff of the metal."""
    d = np.linalg.norm(positions - positions[mi], axis=1)
    return {
        i: float(d[i])
        for i, s in enumerate(symbols)
        if i != mi and s in DONOR_ELEMENTS and d[i] <= cutoff
    }


def _energy(atoms, calc) -> float:
    atoms.calc = calc
    return float(atoms.get_potential_energy())


def _interaction_energy(atoms, calc, mi: int, e_complex: float) -> float | None:
    """Crude metal binding energy: E(complex) − E(apo, frozen) − E(metal). Ranking
    only — no apo relaxation, no solvation, reference-state-naive. Reuses the
    already-relaxed complex energy rather than recomputing it."""
    keep = [i for i in range(len(atoms)) if i != mi]
    try:
        return e_complex - _energy(atoms[keep], calc) - _energy(atoms[[mi]], calc)
    except Exception:
        return None


def relax_site(
    atoms,
    calc,
    metal: str = "Ni",
    cutoff: float = 2.8,
    fmax: float = 0.05,
    steps: int = 200,
    interaction: bool = True,
) -> SiteRelaxation:
    """Relax `atoms` under `calc` and measure how the metal site moved."""
    from ase.optimize import LBFGS

    symbols = atoms.get_chemical_symbols()
    mi = _metal_index(symbols, metal)
    start = atoms.get_positions().copy()
    pre = _shell(start, symbols, mi, cutoff)

    atoms.calc = calc
    LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    energy = float(atoms.get_potential_energy())  # cached at the converged geometry

    pos = atoms.get_positions()
    post = _shell(pos, symbols, mi, cutoff)
    site = [mi, *pre]  # metal + atoms that started in the first shell
    drift = float(np.linalg.norm(pos[site] - start[site], axis=1).max())
    mean_bond = float(np.mean(list(post.values()))) if post else 0.0
    de = _interaction_energy(atoms, calc, mi, energy) if interaction else None

    return SiteRelaxation(len(pre), len(post), drift, mean_bond, de)


class MLIPVerifier:
    """Verifier protocol over an MLIP relaxation. Trusts a design whose metal site
    holds its coordination under the potential; defers when the site collapses or
    the relaxation cannot run. The backbone is pluggable and lazily constructed.
    """

    def __init__(
        self,
        backbone: str = "mace_mp",
        *,
        calculator=None,  # inject a ready ASE calculator (tests, or a shared instance)
        metal_element: str | None = None,  # override; default derives from site.metal
        cutoff: float = 2.8,
        radius: float = 5.0,
        trust_drift: float = 0.5,
        ood_drift: float = 1.5,
        device: str = "cuda",
    ):
        self._backbone = backbone
        self._calc = calculator
        self.metal_element = metal_element
        self.cutoff = cutoff
        self.radius = radius
        self.trust_drift = trust_drift
        self.ood_drift = ood_drift
        self.device = device

    def _metal(self, design: BinderDesign) -> str:
        """ASE element for this design — its site metal, unless overridden."""
        return self.metal_element or element_symbol(design.site.metal)

    @property
    def calc(self):
        if self._calc is None:
            self._calc = make_backbone(self._backbone, self.device)
        return self._calc

    def _cluster(self, design: BinderDesign):
        from ase.io import read

        if not design.source:
            raise ValueError("MLIPVerifier needs design.source (a structure file)")
        atoms = read(design.source)
        pos = atoms.get_positions()
        mi = _metal_index(atoms.get_chemical_symbols(), self._metal(design))
        return atoms[np.linalg.norm(pos - pos[mi], axis=1) <= self.radius]

    def relax(self, design: BinderDesign) -> SiteRelaxation:
        return relax_site(self._cluster(design), self.calc, self._metal(design), self.cutoff)

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            r = self.relax(design)
        except Exception as e:  # divergence / NaN forces / unreadable input ⇒ defer
            return Verdict(0.0, trust=False, ood=True, reason=f"MLIP relaxation failed: {type(e).__name__}")

        held = r.donors_lost == 0
        ood = r.site_drift > self.ood_drift or r.cn_after == 0
        trust = held and r.site_drift <= self.trust_drift and not ood
        # Higher = more stable: decays with site drift, scaled by the fraction of
        # the first shell retained. Interaction energy rides along for ranking.
        score = float(np.exp(-r.site_drift)) * (r.cn_after / max(r.cn_before, 1))
        de = "" if r.interaction_energy is None else f", ΔE_bind {r.interaction_energy:.2f} eV"
        lost = "held" if held else f"lost {r.donors_lost} donor(s)"
        reason = f"site {lost}, drift {r.site_drift:.2f} Å{de}" + (" — defer" if ood else "")
        return Verdict(score, trust=trust, ood=ood, reason=reason)
