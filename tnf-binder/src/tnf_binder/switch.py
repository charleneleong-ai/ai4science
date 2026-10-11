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

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from tnf_binder.pka import run_propka
from tnf_binder.verifiers.pka_switch import pka_window

Site = tuple[str, int]


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
