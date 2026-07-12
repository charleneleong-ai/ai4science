"""Physics-based metal selectivity, via the MLIP binding energy.

The geometry tier can't discriminate divalent metals — Ni/Cu/Co coordination distributions
overlap, so `selectivity_profile` trusts a site for all of them. The MLIP promises a sharper
handle: recompute the metal binding energy ΔE with the site's metal swapped to each competitor,
and trust a design only if its *target* metal binds most favourably.

**That promise is only as good as the backbone's metal physics, and most MLIPs don't have it.**
Metal preference is set by ligand-field stabilisation — d-electron count and spin state (the
Irving–Williams series: Mn<Fe<Co<Ni<Cu>Zn, with Cu²⁺ the peak, for essentially any N/O donor
set). MACE-MP is trained on inorganic crystals, carries no charge or spin state, and empirically
*inverts* the series — it ranks Mn²⁺ strongest and Cu²⁺ fourth. Its metal ordering is therefore
meaningless, and a reward built on it would steer a generator toward nonsense.

So this tier does two things a naive swap-and-relax does not:
  1. **Swaps the spin with the metal** — each divalent 3d ion has a different d-electron count
     (Mn²⁺ d⁵ … Zn²⁺ d¹⁰), so changing only the element compares different electronic systems.
  2. **Gates itself on the Irving–Williams series** (`ranks_irving_williams`) — a backbone that
     cannot reproduce the Cu²⁺ peak on a reference site is refused, and the tier `defer`s rather
     than emit a ranking it cannot justify.

Reuses the MLIP plumbing (`MLIPBase`, `relax_site`) — same backbone and cluster extraction as
the static/dynamics verifiers.
"""

from __future__ import annotations

from dataclasses import dataclass
from weakref import WeakKeyDictionary

import numpy as np

from ..core import BinderDesign, Verdict, element_symbol
from .mlip import MLIPBase, multiplicity, relax_site, single_point

MARGIN_SCALE = 0.3  # eV; fixes the score's sensitivity, independent of the trust cutoff

# The Irving–Williams series, weakest → strongest, peaking at Cu2+ (Zn2+ drops back).
IRVING_WILLIAMS = ("Mn2+", "Fe2+", "Co2+", "Ni2+", "Cu2+", "Zn2+")
# Per-calculator validity, so the probe runs once per backbone. Keyed weakly on the calculator
# itself, never on id(): ids are reused after GC, and a safety gate must not let a freed
# valid backbone hand its verdict to the invalid one that reuses its address.
GATE_CACHE: WeakKeyDictionary = WeakKeyDictionary()


@dataclass
class SelectivityProfile:
    """Per-metal binding energy ΔE (eV) for one design's site. Lower = stronger."""

    target: str
    energies: dict[str, float]

    @property
    def preferred(self) -> str:
        return min(self.energies, key=self.energies.get)

    @property
    def margin(self) -> float:
        """How much more favourably the target binds than the best competitor (eV).
        Positive ⇒ target preferred by that margin; negative ⇒ a competitor wins."""
        competitors = [e for m, e in self.energies.items() if m != self.target]
        return (min(competitors) - self.energies[self.target]) if competitors else 0.0


def swap_metal(atoms, from_el: str, to_el: str, spin: int | None = None):
    """Swap the metal element — *and its spin state*. Each divalent 3d ion has a different
    d-electron count, so swapping the element alone would compare different electronic systems."""
    a = atoms.copy()
    syms = a.get_chemical_symbols()
    syms[syms.index(from_el)] = to_el
    a.set_chemical_symbols(syms)
    if spin is not None:
        a.info["spin"] = spin
    return a


def hexaaqua(metal_el: str):
    """[M(H2O)6]²⁺ — the probe for the validity gate, and not an arbitrary one: octahedral,
    charge +2, closed-shell neutral ligands, and *the very system the Irving–Williams series was
    measured on*. A bare M/N/O cluster would be a set of radicals — meaningless to ask a
    charge/spin-aware backbone to evaluate."""
    from ase import Atoms

    axes = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    symbols, pos = [metal_el], [np.zeros(3)]
    for ax in axes:
        o = ax * 2.10  # M–O
        ref = np.array([0.0, 0.0, 1.0]) if abs(ax[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        perp = np.cross(ax, ref)
        perp /= np.linalg.norm(perp)
        symbols += ["O", "H", "H"]  # H point away from the metal; |O–H| ≈ 0.95 Å
        pos += [o, o + ax * 0.30 + perp * 0.90, o + ax * 0.30 - perp * 0.90]
    return Atoms(symbols=symbols, positions=np.array(pos))


def ranks_irving_williams(calc) -> bool:
    """Can this backbone rank the divalent 3d series at all?

    Probes it on [M(H2O)6]²⁺ and checks it puts **Cu²⁺ at the peak** — the Irving–Williams series,
    the most robust trend in coordination chemistry. A backbone with no charge/spin state (MACE-MP)
    cannot represent the ligand-field stabilisation that drives it, and empirically *inverts* the
    ordering; its metal verdicts are then meaningless and must not be used.

    The probe carries its own state (charge +2, spin from the ion's d-count) — independent of any
    design's cluster charge. E(apo) is the same six waters for every metal, so it cancels from the
    ranking and is never computed. Charge-blind backbones simply ignore the info dict."""
    try:
        cached = GATE_CACHE.get(calc)
        if cached is not None:
            return cached
    except TypeError:  # calculator isn't weak-referenceable ⇒ just re-probe, never mis-cache
        pass

    def remember(ok: bool) -> bool:
        try:
            GATE_CACHE[calc] = ok
        except TypeError:
            pass
        return ok

    energies: dict[str, float] = {}
    for m in IRVING_WILLIAMS:
        mult = multiplicity(m) or 1
        complex_ = hexaaqua(element_symbol(m))
        complex_.info["charge"], complex_.info["spin"] = 2, mult
        ion = complex_[[0]]
        ion.info["charge"], ion.info["spin"] = 2, mult
        try:
            energies[m] = single_point(complex_, calc) - single_point(ion, calc)
        except Exception:  # backbone can't even evaluate the probe ⇒ unusable for ranking
            return remember(False)
    return remember(min(energies, key=energies.get) == "Cu2+")


class MLIPSelectivityVerifier(MLIPBase):
    """Trusts a design whose target metal binds most favourably (by ΔE) among the competitors —
    the discrimination geometry can't make. `defer`s outright if the backbone fails the
    Irving–Williams gate: a metal ranking from a spin-blind potential is not a weak signal, it is
    a meaningless one."""

    def __init__(
        self,
        backbone: str = "mace_mp",
        *,
        metals: tuple[str, ...] = ("Ni2+", "Cu2+", "Co2+"),
        trust_margin: float = 0.2,  # eV the target must beat the best competitor by
        **kw,
    ):
        super().__init__(backbone, **kw)
        self.metals = metals
        self.trust_margin = trust_margin

    def profile(self, design: BinderDesign) -> SelectivityProfile:
        base = self._cluster(design)
        from_el = self._metal(design)
        # always evaluate the design's own metal, even if it isn't in the competitor
        # panel — else SelectivityProfile.margin KeyErrors on the target lookup.
        metals = tuple(dict.fromkeys((design.site.metal, *self.metals)))
        energies = {
            m: relax_site(
                swap_metal(base, from_el, element_symbol(m), spin=multiplicity(m)), self.calc,
                metal=element_symbol(m), interaction=True,
            ).interaction_energy
            for m in metals
        }
        return SelectivityProfile(design.site.metal, energies)

    def verify(self, design: BinderDesign) -> Verdict:
        if not ranks_irving_williams(self.calc):
            return Verdict.defer(
                "backbone cannot rank the divalent 3d series (fails Irving–Williams: Cu2+ is not "
                "the peak) — it has no ligand-field/spin physics, so its metal ordering is meaningless"
            )
        try:
            prof = self.profile(design)
        except Exception as e:  # relaxation / parse failure ⇒ can't judge selectivity
            return Verdict.defer(f"MLIP selectivity failed: {type(e).__name__}")
        if any(e is None for e in prof.energies.values()):
            return Verdict.defer("interaction energy unavailable")

        margin = prof.margin
        trust = prof.preferred == prof.target and margin >= self.trust_margin
        score = float(1.0 / (1.0 + np.exp(-margin / MARGIN_SCALE)))  # 0.5 at margin 0
        reason = f"ΔE favours {prof.preferred} (target {prof.target}, margin {margin:+.2f} eV)"
        return Verdict(score, trust=trust, ood=False, reason=reason)
