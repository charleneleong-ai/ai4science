"""Engineering a pH switch: put histidines somewhere and find out what their pKa actually is.

Three mechanistic rules for *where* a switchable histidine should go were argued for and each
was refuted on real designs — burial lowers pKa rather than raising it (r = -0.37), a
neighbouring cation has the wrong sign (r = +0.11), carboxylate proximity is weak and has
counter-examples (r = +0.21). `verifiers.pka_switch` records that. What is left is to build the
mutant properly and ask a physics solver, which is what this does.

The sequence is always the same and was written out twice: mutate the chosen positions on a
real pose, write the structure, run PROPKA over it, then read back the pKa of each position
that was changed. Both steps have a silent failure mode — a mutation that does not take, and a
residue PROPKA declines to score — so each position reports its own measurement or `None`, and
a caller can tell which positions were answered from which were not.

Single substitutions and combinations are the same operation with one or several positions:
pKa values do not compose (A35H+A38H drops A35 from 6.07 to 5.94, out of the window), so a
combination has to be built and measured as a combination rather than inferred from singles.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypeVar

from tnf_binder.pka import run_propka
from tnf_binder.verifiers.pka_switch import pka_window

Site = tuple[str, int]
# A combination is enumerated over whatever shape the caller carries its positions in: the pKa
# scans use (chain, number), the campaign builder keeps the measured pKa and verdict alongside.
T = TypeVar("T")

ASSAY_PH = (7.4, 6.0)
# Delta-Q is largest for a histidine sitting exactly between the two assay pH values, so what
# costs a design is distance from this midpoint in either direction -- not a low pKa. A pKa of
# 7.6 and one of 5.8 are equally bad bets, which is why the window is symmetric about it.
MIDPOINT = sum(ASSAY_PH) / 2


def protonated_fraction(pka: float, ph: float) -> float:
    """Henderson-Hasselbalch: the fraction of an imidazole carrying a proton at this pH."""
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def delta_q(pka: float) -> float:
    """Charge one histidine gains when the assay moves from pH 7.4 to pH 6.0.

    This is the whole conditional-binding mechanism reduced to one number: at pH 6.0 the
    histidine protonates, the interface gains positive charge it was not designed around, and
    binding should fail. It is a charge, not a free energy -- nothing here claims how much
    affinity a given charge costs, only that more of it at the interface is more disruption.
    """
    high, low = ASSAY_PH
    return protonated_fraction(pka, low) - protonated_fraction(pka, high)



def his_combinations(sites: Sequence[T], max_his: int, cap: int | None = None,
                     ) -> list[tuple[T, ...]]:
    """Every subset of `sites` up to `max_his` positions, smallest first, capped at `cap`.

    Enumerated rather than reasoned about: pKas do not compose, so there is no way to predict
    which subsets are worth building from the singles. Smallest-first ordering matters when a
    cap bites -- each added histidine costs interface confidence (mean binder-to-target ipTM
    0.235 for singles, 0.180 for doubles, 0.158 for triples), so a truncated list should keep
    the cheap combinations and drop the expensive ones.
    """
    combos = [combo for k in range(1, min(max_his, len(sites)) + 1)
              for combo in itertools.combinations(sites, k)]
    return combos if cap is None else combos[:cap]


def his_tag(positions: Sequence[Site]) -> str:
    """The name a mutant is known by: `A35H_A72H`."""
    return "_".join(f"{chain}{number}H" for chain, number in positions)


def position_list(positions: Sequence[Site]) -> str:
    """The positions as a table cell: `A35 A72`."""
    return " ".join(f"{chain}{number}" for chain, number in positions)


@dataclass(frozen=True)
class HisMeasurement:
    """One introduced histidine. `pka` is None when PROPKA declined to score it."""

    chain: str
    resnum: int
    pka: float | None = None

    @property
    def window(self) -> str | None:
        """`ideal` / `usable` / `dead`, or None when there is no measurement to classify."""
        return None if self.pka is None else pka_window(self.pka)

    @property
    def in_window(self) -> bool:
        return self.window in ("ideal", "usable")

    def __str__(self) -> str:
        return (f"{self.chain}{self.resnum}="
                f"{'NA' if self.pka is None else format(self.pka, '.2f')}")


@dataclass(frozen=True)
class HisVariant:
    """A built mutant and what its introduced histidines measured."""

    parent: str
    structure: Path
    measurements: tuple[HisMeasurement, ...]

    @property
    def positions(self) -> tuple[Site, ...]:
        return tuple((m.chain, m.resnum) for m in self.measurements)

    @property
    def name(self) -> str:
        return f"{self.parent}__{his_tag(self.positions)}"

    @property
    def measured(self) -> tuple[HisMeasurement, ...]:
        """Only the histidines PROPKA actually scored."""
        return tuple(m for m in self.measurements if m.pka is not None)

    @property
    def in_window(self) -> int:
        return sum(m.in_window for m in self.measurements)

    @property
    def ideal(self) -> int:
        return sum(m.window == "ideal" for m in self.measurements)


@dataclass(frozen=True)
class SampledHis:
    """One histidine scored against every diffusion sample of its complex, not just one.

    The pKa in the first submission came from a single predicted structure each time. Re-scoring
    against all three Boltz samples moved 8 of 23 histidines out of the window -- R21H alone
    ranges 3.99 to 6.70 where its design structure said 7.05 -- so a single-structure pKa says
    almost nothing on its own and the spread is the measurement.
    """

    chain: str
    resnum: int
    pkas: tuple[float, ...]

    @property
    def delta_q_floor(self) -> float:
        """Delta-Q at the least favourable sample: a lower bound, never an average.

        The worst sample is the one furthest from `MIDPOINT`, which can be the highest pKa as
        easily as the lowest. Taking the minimum pKa instead would flatter any histidine whose
        spread runs above the midpoint.

        Zero when nothing was scored: a histidine PROPKA declined to score contributes no
        measured charge, and must not take the whole ranking down with it.
        """
        return min((delta_q(pka) for pka in self.pkas), default=0.0)

    @property
    def in_window_count(self) -> int:
        """How many of the stored pKas sit in the usable window.

        Named apart from the `in_window_samples` column of the robustness CSV on purpose:
        `switch_profiles` stores only the interval endpoints, so this counts out of 2 where the
        file counts out of `n_samples`. Both answer `robust` identically -- the window is an
        interval, so endpoints inside it put every sample inside it -- but the numbers differ.
        """
        return sum(pka_window(pka) in ("ideal", "usable") for pka in self.pkas)

    @property
    def robust(self) -> bool:
        """True only when every sample agrees the switch is in the usable window."""
        return bool(self.pkas) and self.in_window_count == len(self.pkas)

    @property
    def spread(self) -> float:
        return max(self.pkas) - min(self.pkas) if self.pkas else 0.0


@dataclass(frozen=True)
class SwitchProfile:
    """Every histidine on one design, as the ranking actually saw them."""

    design: str
    histidines: tuple[SampledHis, ...]

    @property
    def delta_q_floor(self) -> float:
        """Charge the whole interface gains at pH 6.0, each histidine at its worst sample.

        Summed rather than averaged: two weak histidines really do perturb an interface more
        than one weak histidine, and the parent with no histidines has to score zero.
        """
        return sum(h.delta_q_floor for h in self.histidines)

    @property
    def robust(self) -> int:
        return sum(h.robust for h in self.histidines)

    @property
    def robust_tag(self) -> str:
        """`2/3` -- how the evidence table reports it."""
        return f"{self.robust}/{len(self.histidines)}"


def switch_profiles(rows: Iterable[dict[str, str]]) -> list[SwitchProfile]:
    """Group robustness rows (one per histidine per design) into one profile per design.

    Reads the long-format CSV `pka_robustness.py` writes, where each row is one histidine with
    its min, max and mean across samples. Only min and max are used: delta-Q is evaluated at the
    extremes because the midpoint-distance rule makes either one the worst case.
    """
    by_design: dict[str, list[SampledHis]] = {}
    for row in rows:
        pkas = tuple(float(row[key]) for key in ("pka_min", "pka_max"))
        by_design.setdefault(row["design"], []).append(
            SampledHis(row.get("chain", "X"), int(row["position"]), pkas))
    return [SwitchProfile(design, tuple(sorted(his, key=lambda h: h.resnum)))
            for design, his in by_design.items()]


def measure_his(parent: Path, positions: Sequence[Site], *, mutator, workdir: Path,
                dest: Path | None = None,
                propka: Callable[[Path, Path], dict[Site, float]] = run_propka,
                ) -> HisVariant | None:
    """Mutate `positions` to histidine on `parent` and measure each one's pKa.

    None when the mutation itself did not take — distinct from a mutant whose pKa could not be
    predicted, which comes back as a variant whose measurements are None. The caller needs both
    cases separable: the first is a broken build, the second is a question PROPKA declined.

    `dest` is where the mutant is kept when it outlives the measurement -- a combination that
    will be refolded later, say -- and defaults to the scratch directory PROPKA works in.

    `propka` is injectable so the sequence can be tested without the solver installed.
    """
    pose = mutator.load(parent)
    if not mutator.mutate(pose, list(positions)):
        return None
    mutant = (dest or workdir) / f"{parent.stem}__{his_tag(positions)}.pdb"
    mutator.write(pose, mutant)
    pkas = propka(mutant, workdir)
    return HisVariant(parent.stem, mutant,
                      tuple(HisMeasurement(chain, number, pkas.get((chain, number)))
                            for chain, number in positions))
