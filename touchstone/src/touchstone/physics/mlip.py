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


@dataclass
class SiteDynamics:
    """Metal-site coordination across a short MLIP molecular-dynamics run."""

    cn_initial: int
    retention: float  # fraction of sampled frames that kept the full first shell


def md_site(
    atoms,
    calc,
    metal: str = "Ni",
    cutoff: float = 2.8,
    temperature: float = 300.0,
    steps: int = 500,
    timestep: float = 1.0,
    friction: float = 0.02,
    sample_every: int = 10,
) -> SiteDynamics:
    """Run short NVT (Langevin) MD and measure how often the first shell survives."""
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    symbols = atoms.get_chemical_symbols()
    mi = _metal_index(symbols, metal)
    cn0 = len(_shell(atoms.get_positions(), symbols, mi, cutoff))

    atoms.calc = calc
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
    dyn = Langevin(atoms, timestep * units.fs, temperature_K=temperature, friction=friction)
    kept: list[bool] = []

    def sample():
        kept.append(len(_shell(atoms.get_positions(), symbols, mi, cutoff)) >= cn0)

    dyn.attach(sample, interval=sample_every)
    dyn.run(steps)
    return SiteDynamics(cn0, sum(kept) / len(kept) if kept else 0.0)


class _MLIPBase:
    """Shared plumbing for MLIP verifiers: a lazily-built, pluggable backbone and
    the metal-centred cluster pulled from a design's structure file.

    Pass one `make_backbone(...)` instance as `calculator=` to both the static and
    dynamics verifiers to share a single in-memory model across them.
    """

    def __init__(
        self,
        backbone: str = "mace_mp",
        *,
        calculator=None,  # inject a ready ASE calculator (tests, or a shared instance)
        metal_element: str | None = None,  # override; default derives from site.metal
        cutoff: float = 2.8,
        radius: float = 5.0,
        device: str = "cuda",
    ):
        self._backbone = backbone
        self._calc = calculator
        self.metal_element = metal_element
        self.cutoff = cutoff
        self.radius = radius
        self.device = device

    @property
    def calc(self):
        if self._calc is None:
            self._calc = make_backbone(self._backbone, self.device)
        return self._calc

    def _metal(self, design: BinderDesign) -> str:
        """ASE element for this design — its site metal, unless overridden."""
        return self.metal_element or element_symbol(design.site.metal)

    def _cluster(self, design: BinderDesign):
        from ase.io import read

        if not design.source:
            raise ValueError("MLIP verifier needs design.source (a structure file)")
        atoms = read(design.source)
        pos = atoms.get_positions()
        mi = _metal_index(atoms.get_chemical_symbols(), self._metal(design))
        return atoms[np.linalg.norm(pos - pos[mi], axis=1) <= self.radius]


class MLIPVerifier(_MLIPBase):
    """Trusts a design whose metal site holds its coordination under an MLIP
    relaxation; defers when the site collapses or the relaxation cannot run."""

    def __init__(self, backbone: str = "mace_mp", *, trust_drift: float = 0.5, ood_drift: float = 1.5, **kw):
        super().__init__(backbone, **kw)
        self.trust_drift = trust_drift
        self.ood_drift = ood_drift

    def relax(self, design: BinderDesign) -> SiteRelaxation:
        return relax_site(self._cluster(design), self.calc, self._metal(design), self.cutoff)

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            r = self.relax(design)
        except Exception as e:  # divergence / NaN forces / unreadable input ⇒ defer
            return Verdict.defer(f"MLIP relaxation failed: {type(e).__name__}")

        held = r.donors_lost == 0
        # Higher = more stable: decays with site drift, scaled by the fraction of
        # the first shell retained. Clamp to [0,1] — a donor migrating into the shell
        # can push cn_after/cn_before above 1. Interaction energy rides along for ranking.
        score = min(1.0, float(np.exp(-r.site_drift)) * (r.cn_after / max(r.cn_before, 1)))
        de = "" if r.interaction_energy is None else f", ΔE_bind {r.interaction_energy:.2f} eV"
        lost = "held" if held else f"lost {r.donors_lost} donor(s)"
        reason = f"site {lost}, drift {r.site_drift:.2f} Å{de}"
        if r.site_drift > self.ood_drift or r.cn_after == 0:
            return Verdict.defer(reason, score=score)
        return Verdict(score, trust=held and r.site_drift <= self.trust_drift, ood=False, reason=reason)


class MLIPDynamicsVerifier(_MLIPBase):
    """Trusts a design whose coordination survives short MLIP molecular dynamics —
    a thermal-stability check. The dynamic counterpart to MLIPVerifier's static
    relax, and an independent second method to the xtb cluster-MD tier."""

    def __init__(
        self,
        backbone: str = "mace_mp",
        *,
        temperature: float = 300.0,
        trust_retention: float = 0.8,
        ood_retention: float = 0.5,
        steps: int = 500,
        **kw,
    ):
        super().__init__(backbone, **kw)
        self.temperature = temperature
        self.trust_retention = trust_retention
        self.ood_retention = ood_retention
        self.steps = steps

    def dynamics(self, design: BinderDesign) -> SiteDynamics:
        return md_site(
            self._cluster(design), self.calc, self._metal(design),
            self.cutoff, self.temperature, self.steps,
        )

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            d = self.dynamics(design)
        except Exception as e:  # blow-up / unreadable input ⇒ defer
            return Verdict.defer(f"MLIP MD failed: {type(e).__name__}")

        if d.cn_initial == 0:  # a naked metal trivially "retains" 0 donors — don't trust it
            return Verdict.defer("no coordinating atoms to track")
        # retention is higher-is-better, so the bound direction inverts vs the
        # strain/drift verifiers: trust above trust_retention, defer below ood_retention.
        reason = f"shell survived {d.retention:.0%} of {self.temperature:.0f} K MD"
        if d.retention < self.ood_retention:
            return Verdict.defer(reason, score=d.retention)
        return Verdict(d.retention, trust=d.retention >= self.trust_retention, ood=False, reason=reason)
