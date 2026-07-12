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

import atexit
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core import BinderDesign, Verdict, element_symbol
from ..geometry.parse import DONOR_ELEMENTS

H_CACHE: dict[tuple[str, float], str] = {}  # (source, pH) → protonated path; dedups across verifiers/batches
H_DIR: str | None = None  # one temp dir for all protonated copies, removed at process exit


def h_dir() -> str:
    global H_DIR
    if H_DIR is None:
        H_DIR = tempfile.mkdtemp(prefix="touchstone_h_")
        atexit.register(shutil.rmtree, H_DIR, ignore_errors=True)
    return H_DIR


def protonate(structure: str | Path, pH: float = 7.4) -> str | None:
    """Add hydrogens to a structure for the MLIP tier and return the path to a
    protonated copy (or None if OpenBabel isn't installed — the caller then proceeds
    unprotonated). MACE needs explicit H: on a bare BoltzGen backbone it sees
    under-coordinated His/backbone heavy atoms and the metal wanders out of its site.
    OpenBabel at `pH` leaves metal-coordinating donors (N/S bonded to the cation)
    deprotonated, so the donating lone pairs are preserved. Heavy-atom coordinates are
    untouched — only H are added — so the geometry/bond-valence verdicts are unchanged.

    Cached per (structure, pH): the static and dynamics verifiers share one protonation,
    and copies live in a single temp dir cleaned up at process exit."""
    try:
        from openbabel import pybel
    except ImportError:  # optional MLIP-tier dep; degrade to no protonation
        return None

    structure = Path(structure)
    key = (str(structure), pH)
    if key not in H_CACHE:
        mol = next(pybel.readfile("pdb", normalize(structure)))
        mol.OBMol.AddHydrogens(False, True, pH)  # (polar-only=False, correct-for-pH=True, pH)
        out = str(Path(h_dir()) / f"{len(H_CACHE)}_{structure.stem}_H.pdb")
        mol.write("pdb", out, overwrite=True)
        H_CACHE[key] = out
    return H_CACHE[key]


def normalize(structure: Path) -> str:
    """Rewrite to a clean PDB via gemmi so OpenBabel/ASE can type the metal — BoltzGen
    mmCIFs encode the cation as element `*`/unknown (→ downstream KeyError); gemmi reads it.
    Falls back to the raw path if gemmi is absent (fine for PDBs; CIF metals may be lost)."""
    try:
        import gemmi
    except ImportError:
        return str(structure)
    st = gemmi.read_structure(str(structure))
    out = str(Path(h_dir()) / f"{structure.stem}_norm.pdb")
    st.write_pdb(out)
    return out


def make_backbone(name: str, device: str = "cuda"):
    """ASE calculator for a backbone. Heavy imports are deferred to call time."""
    if name == "mace_mp":
        from mace.calculators import mace_mp

        return mace_mp(model="medium", device=device, default_dtype="float64")
    if name == "uma":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        unit = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)
        return FAIRChemCalculator(unit, task_name="omol")
    if name == "orbmol":  # non-equivariant (no e3nn) ⇒ co-resides with MACE in one env
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator

        orbff, adapter = pretrained.orbmol_v2(device=device, precision="float64")
        return ORBCalculator(orbff, atoms_adapter=adapter, device=device)
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
    converged: bool = True  # optimiser reached fmax within the step budget
    max_force: float = 0.0  # final |F|max over the free atoms (eV/Å) — the relaxation-health signal

    @property
    def donors_lost(self) -> int:
        return max(0, self.cn_before - self.cn_after)


def metal_index(symbols: list[str], metal: str) -> int:
    if metal not in symbols:
        raise ValueError(f"no {metal} atom in cluster")
    return symbols.index(metal)


def shell(positions: np.ndarray, symbols: list[str], mi: int, cutoff: float) -> dict[int, float]:
    """{donor index: distance} for N/O/S within cutoff of the metal."""
    d = np.linalg.norm(positions - positions[mi], axis=1)
    return {
        i: float(d[i])
        for i, s in enumerate(symbols)
        if i != mi and s in DONOR_ELEMENTS and d[i] <= cutoff
    }


def potential_energy(atoms, calc) -> float:
    atoms.calc = calc
    return float(atoms.get_potential_energy())


def freeze_scaffold(atoms, free: set[int]) -> None:
    """Position-restrain every atom outside `free` (the metal + first shell). A
    metal-centred cluster is a cut-out of a protein: without the backbone holding the
    donors, free relaxation lets the whole fragment disperse (donors drift many Å,
    energies blow up). Freezing the scaffold lets only the coordination polyhedron
    relax against a fixed backbone — the standard frozen-boundary cluster treatment."""
    from ase.constraints import FixAtoms

    frozen = [i for i in range(len(atoms)) if i not in free]
    if frozen:
        atoms.set_constraint(FixAtoms(indices=frozen))


def single_point(atoms, calc) -> float:
    """Constraint-free single-point energy of a freshly-sliced cluster. Callers pass
    `atoms[...]` slices (already independent copies); a slice carries stale FixAtoms
    indices, so clear them in place — irrelevant to a single point anyway."""
    atoms.set_constraint()
    return potential_energy(atoms, calc)


def interaction_energy(atoms, calc, mi: int, e_complex: float) -> float | None:
    """Crude metal binding energy: E(complex) − E(apo, frozen) − E(metal). Ranking
    only — no apo relaxation, no solvation, reference-state-naive. Reuses the
    already-relaxed complex energy rather than recomputing it."""
    keep = [i for i in range(len(atoms)) if i != mi]
    try:
        return e_complex - single_point(atoms[keep], calc) - single_point(atoms[[mi]], calc)
    except (RuntimeError, ValueError, FloatingPointError):
        return None


def relax_site(
    atoms,
    calc,
    metal: str = "Ni",
    cutoff: float = 2.8,
    fmax: float = 0.05,
    steps: int = 200,
    interaction: bool = True,
    restrain: bool = True,
) -> SiteRelaxation:
    """Relax `atoms` under `calc` and measure how the metal site moved. With
    `restrain` (default), the backbone scaffold is frozen so only the metal + first
    shell relax — without it a cut-out cluster disperses unphysically."""
    from ase.optimize import LBFGS

    symbols = atoms.get_chemical_symbols()
    mi = metal_index(symbols, metal)
    start = atoms.get_positions().copy()
    pre = shell(start, symbols, mi, cutoff)

    if restrain:
        freeze_scaffold(atoms, {mi, *pre})
    atoms.calc = calc
    converged = bool(LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps))
    energy = float(atoms.get_potential_energy())  # cached at the converged geometry
    if not np.isfinite(energy):  # diverged into non-finite territory ⇒ unjudgeable
        raise FloatingPointError("non-finite relaxation energy")
    # final |F|max over the free atoms (FixAtoms zeroes the frozen ones) — the same
    # basis the optimiser's convergence test uses, kept as a graded relaxation-health metric.
    max_force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())

    pos = atoms.get_positions()
    post = shell(pos, symbols, mi, cutoff)
    site = [mi, *pre]  # metal + atoms that started in the first shell
    drift = float(np.linalg.norm(pos[site] - start[site], axis=1).max())
    mean_bond = float(np.mean(list(post.values()))) if post else 0.0
    de = interaction_energy(atoms, calc, mi, energy) if interaction else None

    return SiteRelaxation(len(pre), len(post), drift, mean_bond, de, converged, max_force)


@dataclass
class SitePreorganization:
    """How the coordination shell reorganizes when the metal is removed (apo) and relaxed."""

    donor_drift: float  # max first-shell donor displacement from the holo positions (Å)
    donors: int  # first-shell donor count in the holo


def relax_apo(
    atoms, calc, metal: str = "Ni", cutoff: float = 2.8, fmax: float = 0.05, steps: int = 200,
) -> SitePreorganization:
    """Remove the metal, relax the freed first-shell donors against the frozen backbone, and
    return the max donor drift from the holo positions (the metal-off reorganization)."""
    from ase.optimize import LBFGS

    symbols = atoms.get_chemical_symbols()
    mi = metal_index(symbols, metal)
    start = atoms.get_positions().copy()
    pre = list(shell(start, symbols, mi, cutoff))  # holo first-shell donor indices
    if not pre:
        raise ValueError("no first-shell donors to track")
    holo = start[pre]  # holo donor positions
    apo = atoms[[i for i in range(len(atoms)) if i != mi]]  # independent copy without the metal
    donors = [i - (i > mi) for i in pre]  # reindex donors after dropping the metal
    freeze_scaffold(apo, set(donors))  # only the donors relax, against the frozen backbone
    apo.calc = calc
    LBFGS(apo, logfile=None).run(fmax=fmax, steps=steps)
    drift = float(np.linalg.norm(apo.get_positions()[donors] - holo, axis=1).max())
    return SitePreorganization(drift, len(pre))


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
    restrain: bool = True,
) -> SiteDynamics:
    """Run short NVT (Langevin) MD and measure how often the first shell survives.
    With `restrain` (default), the backbone scaffold is frozen so the cluster cannot
    disperse — the same frozen-boundary treatment as `relax_site`."""
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    symbols = atoms.get_chemical_symbols()
    mi = metal_index(symbols, metal)
    shell0 = shell(atoms.get_positions(), symbols, mi, cutoff)
    cn0 = len(shell0)

    if restrain:
        freeze_scaffold(atoms, {mi, *shell0})
    atoms.calc = calc
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
    dyn = Langevin(atoms, timestep * units.fs, temperature_K=temperature, friction=friction)
    kept: list[bool] = []

    def sample():
        kept.append(len(shell(atoms.get_positions(), symbols, mi, cutoff)) >= cn0)

    dyn.attach(sample, interval=sample_every)
    dyn.run(steps)
    return SiteDynamics(cn0, sum(kept) / len(kept) if kept else 0.0)


class MLIPBase:
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
        restrain: bool = True,  # freeze the backbone scaffold during relax/MD
        protonate: bool = True,  # add H first (MACE needs it); no-op if OpenBabel absent
        charge: int | None = None,  # total cluster charge — required by charge-aware backbones (OrbMol)
        spin: int | None = None,  # spin multiplicity — required by charge-aware backbones (OrbMol)
    ):
        # OrbMol reads total charge + spin off atoms.info and cannot run without them.
        # Fail loud at construction, not deep in the orb adapter — and total cluster
        # charge can't be honestly derived (it depends on donor protonation), so it must
        # be supplied. Other backbones ignore atoms.info, so the fields stay optional.
        if backbone == "orbmol" and (charge is None or spin is None):
            raise ValueError("orbmol backbone requires charge and spin (total cluster charge + multiplicity)")
        self._backbone = backbone
        self._calc = calculator
        self.metal_element = metal_element
        self.cutoff = cutoff
        self.radius = radius
        self.device = device
        self.restrain = restrain
        self.protonate = protonate
        self.charge = charge
        self.spin = spin

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
        # add H first (MACE needs explicit H or the metal wanders); falls back to the raw
        # structure if OpenBabel is absent. Heavy-atom coords are untouched, so the parsed
        # site (geometry/bond-valence) is unaffected — only this MLIP cluster sees the H.
        # normalize either way (protonate does it internally; the no-protonation fallback
        # needs it too, else a BoltzGen CIF's untyped metal breaks the ASE read below)
        source = (self.protonate and protonate(design.source)) or normalize(design.source)
        atoms = read(source)
        # PDBs (incl. BoltzGen's) carry a CRYST1 1 Å cell, so ase.io.read marks the
        # structure periodic. Treated as a crystal, an MLIP builds a vast periodic
        # neighbour list (millions of image edges ⇒ OOM). We want a finite cluster.
        atoms.set_pbc(False)
        atoms.set_cell(None)
        pos = atoms.get_positions()
        mi = metal_index(atoms.get_chemical_symbols(), self._metal(design))
        cluster = atoms[np.linalg.norm(pos - pos[mi], axis=1) <= self.radius]
        # Standard ASE per-structure metadata; consumed only by charge/spin-aware
        # backbones (OrbMol), ignored by the rest. Written only when set — absence, not
        # a spurious default, so a MACE/UMA cluster is byte-identical to before.
        if self.charge is not None:
            cluster.info["charge"] = self.charge
        if self.spin is not None:
            cluster.info["spin"] = self.spin
        return cluster


class MLIPVerifier(MLIPBase):
    """Trusts a design whose metal site holds its coordination under an MLIP
    relaxation; defers when the site collapses or the relaxation cannot run."""

    def __init__(
        self, backbone: str = "mace_mp", *, trust_drift: float = 0.5, ood_drift: float = 1.5, steps: int = 200, **kw
    ):
        super().__init__(backbone, **kw)
        self.trust_drift = trust_drift
        self.ood_drift = ood_drift
        self.steps = steps

    def relax(self, design: BinderDesign) -> SiteRelaxation:
        return relax_site(
            self._cluster(design), self.calc, self._metal(design), self.cutoff,
            steps=self.steps, restrain=self.restrain,
        )

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
        metrics = {
            "drift_angstrom": round(r.site_drift, 3), "cn_before": r.cn_before, "cn_after": r.cn_after,
            "max_force_ev_ang": round(r.max_force, 3), "converged": r.converged,
        }
        if r.interaction_energy is not None:
            metrics["interaction_energy_ev"] = round(r.interaction_energy, 3)
        # Relaxation-health gate: a relaxation that never reached fmax leaves drift/CN
        # read off a non-equilibrium geometry — the surrogate's own "I didn't settle"
        # signal, the materials analog of a low confidence score. Defer before drift/CN.
        if not r.converged:
            return Verdict.defer(
                f"relaxation did not settle: |F|max {r.max_force:.2f} eV/Å after {self.steps} steps",
                score=score, metrics=metrics,
            )

        de = "" if r.interaction_energy is None else f", ΔE_bind {r.interaction_energy:.2f} eV"
        lost = "held" if held else f"lost {r.donors_lost} donor(s)"
        reason = f"site {lost}, drift {r.site_drift:.2f} Å{de}"
        if r.site_drift > self.ood_drift or r.cn_after == 0:
            return Verdict.defer(reason, score=score, metrics=metrics)
        return Verdict(score, trust=held and r.site_drift <= self.trust_drift, ood=False, reason=reason, metrics=metrics)


class MLIPDynamicsVerifier(MLIPBase):
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
            self.cutoff, self.temperature, self.steps, restrain=self.restrain,
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
        metrics = {"retention": round(d.retention, 3), "cn_initial": d.cn_initial, "temperature_k": self.temperature}
        if d.retention < self.ood_retention:
            return Verdict.defer(reason, score=d.retention, metrics=metrics)
        return Verdict(d.retention, trust=d.retention >= self.trust_retention, ood=False, reason=reason, metrics=metrics)
