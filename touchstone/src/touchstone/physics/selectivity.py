"""Physics-based metal selectivity, via the MLIP binding energy.

The geometry tier can't discriminate divalent metals — Ni/Cu/Co coordination
distributions overlap, so `selectivity_profile` trusts a site for all of them. The
MLIP gives a sharper handle: recompute the metal binding energy ΔE with the site's
metal swapped to each competitor, and trust a design only if its *target* metal
binds most favourably. Differential ΔE is exactly the donor-identity/HSAB signal
geometry lacks.

Reuses the MLIP plumbing (`_MLIPBase`) and `relax_site` — the same backbone and
cluster extraction as the static/dynamics verifiers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core import BinderDesign, Verdict, element_symbol, reraise_if_bug
from .mlip import _MLIPBase, relax_site

_MARGIN_SCALE = 0.3  # eV; fixes the score's sensitivity, independent of the trust cutoff


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


def _swap_metal(atoms, from_el: str, to_el: str):
    a = atoms.copy()
    syms = a.get_chemical_symbols()
    syms[syms.index(from_el)] = to_el
    a.set_chemical_symbols(syms)
    return a


class MLIPSelectivityVerifier(_MLIPBase):
    """Trusts a design whose target metal binds most favourably (by ΔE) among the
    competitors — the discrimination geometry can't make."""

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
                _swap_metal(base, from_el, element_symbol(m)), self.calc,
                metal=element_symbol(m), interaction=True,
            ).interaction_energy
            for m in metals
        }
        return SelectivityProfile(design.site.metal, energies)

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            prof = self.profile(design)
        except Exception as e:  # relaxation / parse failure ⇒ can't judge selectivity
            reraise_if_bug(e)
            return Verdict.defer(f"MLIP selectivity failed: {e}")
        if any(e is None for e in prof.energies.values()):
            return Verdict.defer("interaction energy unavailable")

        margin = prof.margin
        trust = prof.preferred == prof.target and margin >= self.trust_margin
        score = float(1.0 / (1.0 + np.exp(-margin / _MARGIN_SCALE)))  # 0.5 at margin 0
        reason = f"ΔE favours {prof.preferred} (target {prof.target}, margin {margin:+.2f} eV)"
        return Verdict(score, trust=trust, ood=False, reason=reason)
