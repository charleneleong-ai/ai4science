"""Judge the pH switch on predicted pKa, because reasoning about it failed three times.

Objective 3 wants binding at pH 7.4 and release at pH 6.0, which requires a histidine whose
pKa lies BETWEEN them: neutral at 7.4 so it binds, cationic at 6.0 so it lets go.

Three mechanisms were reasoned out and each was refuted by PROPKA on real designs:

  buried His pays a desolvation penalty       burial LOWERS pKa (r = -0.37); buried interface
                                              positions came back at 3.1-3.6, never firing
  His beside a cation is repelled when        cation count correlates with the WRONG sign
  protonated                                  (r = +0.11)
  carboxylate proximity raises pKa            r = +0.21, with counter-examples at net +6
                                              carboxylates landing at pKa 3.6

No structural descriptor predicted it. Positions therefore come from a per-position scan,
never from a rule -- which is why this verifier consumes measured pKa rather than geometry.
"""

from __future__ import annotations

from dataclasses import dataclass

from touchstone.core import Verdict

from tnf_binder.design import TnfDesign

IDEAL_WINDOW = (6.3, 7.2)
USABLE_WINDOW = (6.0, 7.6)
MIN_INTERFACE_DSASA = 10.0


@dataclass(frozen=True)
class HistidineSite:
    """One histidine: where it sits, how buried it is, and its predicted pKa."""

    position: int
    pka: float
    interface_dsasa: float

    @property
    def at_interface(self) -> bool:
        return self.interface_dsasa >= MIN_INTERFACE_DSASA

    @property
    def window(self) -> str:
        """ideal | usable | dead — where the pKa sits relative to the assay."""
        if IDEAL_WINDOW[0] <= self.pka <= IDEAL_WINDOW[1]:
            return "ideal"
        if USABLE_WINDOW[0] <= self.pka <= USABLE_WINDOW[1]:
            return "usable"
        return "dead"

    @property
    def switches(self) -> bool:
        """A histidine only switches if it is in the window AND at the interface."""
        return self.at_interface and self.window in ("ideal", "usable")


class PkaSwitchVerifier:
    """Trusts a design carrying at least `min_switches` switchable histidines.

    Several in-window histidines give a steeper, cooperative transition than one. Measured
    shifts on combining them were small (<=0.21 pKa units), but combinations are still scored
    as measured rather than assumed to compose.

    Note this predicts protonation state, not binding: no folding model represents pH, and the
    assay is not simulated.
    """

    def __init__(self, *, min_switches: int = 1, prefer_ideal: bool = True) -> None:
        self.min_switches = min_switches
        self.prefer_ideal = prefer_ideal

    def verify(self, design: TnfDesign) -> Verdict:
        sites = design.metrics.get("histidine_sites")
        if sites is None:
            return Verdict.defer("no pKa prediction available")
        if not sites:
            return Verdict(0.0, trust=False, ood=False,
                           reason="no histidines at the interface", metrics={"switches": 0})

        ideal = [s for s in sites if s.switches and s.window == "ideal"]
        usable = [s for s in sites if s.switches and s.window == "usable"]
        dead = [s for s in sites if not s.switches]
        score = len(ideal) + 0.4 * len(usable) if self.prefer_ideal else len(ideal) + len(usable)
        metrics = {
            "ideal": len(ideal), "usable": len(usable), "dead": len(dead),
            "switches": len(ideal) + len(usable),
            "pkas": [f"{s.position}:{s.pka:.2f}" for s in sites],
        }
        trust = (len(ideal) + len(usable)) >= self.min_switches
        detail = ", ".join(f"H{s.position} pKa {s.pka:.2f} ({s.window})" for s in sites[:4])
        return Verdict(round(score, 3), trust=trust, ood=False, metrics=metrics,
                       reason=f"{len(ideal)} ideal + {len(usable)} usable switches — {detail}")
