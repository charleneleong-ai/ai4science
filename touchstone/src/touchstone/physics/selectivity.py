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

So this tier does three things a naive swap-and-relax does not:
  1. **Swaps charge and spin with the metal** (`swap_metal`) — the d-count depends on the ion, not
     the element (Fe²⁺ d⁶ vs Fe³⁺ d⁵), so changing only the element compares different electronic
     systems and leaves the cluster's total charge behind.
  2. **Refuses ions whose spin state it cannot state** — `multiplicity` returns None rather than a
     default. A guessed spin doesn't fail loudly; it silently moves the ligand-field energy, which
     is the number this tier reports.
  3. **Gates itself on the whole Irving–Williams series** (`ranks_irving_williams`) — not just the
     Cu²⁺ peak, which is one bit and easy to hit by accident. A backbone that cannot reproduce the
     series is refused and the tier `defer`s rather than emit a ranking it cannot justify.

**Known gap:** the gate's probe is all-oxygen (hexaaqua), while selectivity in real designs is
largely a donor-identity (HSAB) effect — soft Cys-S/Met-S vs hard N/O. A backbone that passes has
been shown to rank metals on *hard-donor* sites only. See `ranks_irving_williams`.

Reuses the MLIP plumbing (`MLIPBase`, `relax_site`) — same backbone and cluster extraction as
the static/dynamics verifiers.
"""

from __future__ import annotations

from dataclasses import dataclass
from weakref import WeakKeyDictionary

import numpy as np

from ..core import BinderDesign, Verdict, element_symbol, oxidation_state
from .mlip import MLIPBase, multiplicity, relax_site, single_point

MARGIN_SCALE = 0.3  # eV; fixes the score's sensitivity, independent of the trust cutoff

# The Irving–Williams rising limb, weakest → strongest, peaking at Cu2+. Zn2+ is checked separately:
# it must drop back *below* Cu2+ (d10 — no ligand-field stabilisation left), so it isn't part of the
# monotonic run.
IRVING_WILLIAMS = ("Mn2+", "Fe2+", "Co2+", "Ni2+", "Cu2+")
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


def swap_metal(atoms, from_el: str, to_ion: str):
    """Swap the metal — element, oxidation state, charge *and* spin, which travel together.

    Swapping the element alone compares physically different electronic systems (Mn²⁺ d⁵ against
    Zn²⁺ d¹⁰ at the same multiplicity). Swapping element+spin but not charge breaks any panel that
    mixes oxidation states (Cu²⁺ vs Cu⁺, Fe²⁺ vs Fe³⁺): the cluster's total charge must move with
    the ion, or the apo leg silently absorbs the difference."""
    a = atoms.copy()
    syms = a.get_chemical_symbols()
    syms[syms.index(from_el)] = element_symbol(to_ion)
    a.set_chemical_symbols(syms)

    old = a.info.get("ion")
    a.info["ion"] = to_ion
    if (mult := multiplicity(to_ion)) is not None:
        a.info["spin"] = mult
    if a.info.get("charge") is not None and old:
        a.info["charge"] += oxidation_state(to_ion) - oxidation_state(old)
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

    Probes it on [M(H2O)6]²⁺ and requires the **whole Irving–Williams series**: stability rising
    monotonically Mn²⁺ < Fe²⁺ < Co²⁺ < Ni²⁺ < Cu²⁺, then dropping back at Zn²⁺ (d¹⁰ — no
    ligand-field stabilisation left). This is the most robust trend in coordination chemistry, and it
    holds for essentially any N/O donor set.

    Requiring the *series*, not just the Cu²⁺ peak, is deliberate: a single argmin check is one bit,
    and a backbone that puts Cu on top while ranking Mn²⁺ second-strongest would sail through it
    while still being physically meaningless. A backbone with no charge/spin state (MACE-MP) cannot
    represent the ligand-field stabilisation that drives the series at all.

    **Scope — read before trusting a pass.** The probe is all-oxygen (aqua). Selectivity in real
    designs is largely a *donor-identity* (HSAB) effect — soft Cys-S/Met-S vs hard N/O — so passing
    this gate says a backbone can rank metals on hard-donor sites. It does **not** establish that it
    can do so on the thiolate/imidazole sites the tier actually judges, and it says nothing about
    metals outside the divalent 3d series (Pd²⁺/Pt²⁺/Au³⁺ — the recovery targets). Closing that gap
    needs a soft-donor probe with reference data we don't currently have. See
    docs/experiments/2026-07-13-mlip-cannot-rank-metals.md.

    Each probe carries its own state (charge +2, spin from the ion's d-count) — independent of any
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
    for m in (*IRVING_WILLIAMS, "Zn2+"):
        mult = multiplicity(m)
        if mult is None:  # our own table must cover the probe — never guess a spin state
            return remember(False)
        complex_ = hexaaqua(element_symbol(m))
        complex_.info["charge"], complex_.info["spin"] = 2, mult
        ion = complex_[[0]]
        ion.info["charge"], ion.info["spin"] = 2, mult
        try:
            energies[m] = single_point(complex_, calc) - single_point(ion, calc)
        except Exception:  # backbone can't even evaluate the probe ⇒ unusable for ranking
            return remember(False)

    # binding energies: lower = more stable. The rising limb must strictly strengthen Mn→Cu, and
    # Zn2+ must fall back below Cu2+ — the peak alone is one bit and far too easy to satisfy.
    limb = [energies[m] for m in IRVING_WILLIAMS]
    rises = all(a > b for a, b in zip(limb, limb[1:]))
    zn_drops = energies["Zn2+"] > energies["Cu2+"]
    return remember(rises and zn_drops)


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

    def panel(self, design: BinderDesign) -> tuple[str, ...]:
        # always evaluate the design's own metal, even if it isn't in the competitor
        # panel — else SelectivityProfile.margin KeyErrors on the target lookup.
        return tuple(dict.fromkeys((design.site.metal, *self.metals)))

    def profile(self, design: BinderDesign) -> SelectivityProfile:
        base = self._cluster(design)
        from_el = self._metal(design)
        energies = {
            m: relax_site(
                swap_metal(base, from_el, m), self.calc, metal=element_symbol(m), interaction=True,
            ).interaction_energy
            for m in self.panel(design)
        }
        return SelectivityProfile(design.site.metal, energies)

    def verify(self, design: BinderDesign) -> Verdict:
        if unknown := [m for m in self.panel(design) if multiplicity(m) is None]:
            # no tabulated d-count ⇒ no spin state. A default (singlet) wouldn't fail loudly, it
            # would just move the ligand-field energy — which is the number this tier reports.
            return Verdict.defer(
                f"no tabulated spin state for {', '.join(unknown)} — refusing to guess: the spin "
                "state *is* the ligand-field physics that decides metal preference"
            )
        if not ranks_irving_williams(self.calc):
            return Verdict.defer(
                "backbone cannot rank the divalent 3d series (fails Irving–Williams: stability must "
                "rise Mn<Fe<Co<Ni<Cu and drop at Zn) — it has no ligand-field/spin physics, so its "
                "metal ordering is meaningless"
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
