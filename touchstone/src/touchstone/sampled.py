"""A measurement taken on every diffusion sample, and whether the verdict survives the choice.

A structure predictor emits several samples per design. Score one of them and you have a
number; score all of them and you find out whether that number was a property of the design or
of the sample you happened to open. The difference is not academic -- it has decided real
verdicts three times in one campaign:

  - pKa, scored on the design structure, put 8 of 23 histidines in the switching window that
    leave it on a sibling sample. One histidine reads 3.99 and 6.70 on two structures the
    folding model scores identically.
  - interface confidence, taken best-of-n, overstates the mean by up to 0.35 within one design
    and reorders the ranking.
  - steric clashes and buried area, scored on `model_0`, flip the clash-free verdict for 10 of
    16 candidates, with one design reading 0 clashes over 658 A^2 on one sample and 9 clashes
    over 1203 A^2 on another.

Each was found separately and fixed separately, twice, before the third was noticed. This is
the shared part: wrap the per-sample values, and make a gate say whether it would survive being
applied to a different sample.

`Sampled.worst` rather than the mean is the default for a gate, because a design that fails on
any sample has an interpretation the mean does not -- the mean of a clash count over structures
the model cannot distinguish is not a quantity anything measures.
"""

from __future__ import annotations

import statistics as st
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Sampled:
    """One quantity measured on each sample of a design."""

    name: str
    values: tuple[float, ...]

    @classmethod
    def of(cls, name: str, values: Iterable[float]) -> "Sampled":
        return cls(name, tuple(values))

    def __bool__(self) -> bool:
        return bool(self.values)

    @property
    def mean(self) -> float:
        return st.fmean(self.values) if self.values else 0.0

    @property
    def spread(self) -> float:
        return max(self.values) - min(self.values) if self.values else 0.0

    def worst(self, *, higher_is_worse: bool = True) -> float:
        """The least favourable sample, which is what a conservative gate should read."""
        if not self.values:
            return 0.0
        return max(self.values) if higher_is_worse else min(self.values)

    def agreement(self, gate: Callable[[float], bool]) -> float:
        """Fraction of samples on which the gate says keep. 1.0 or 0.0 means they all agree."""
        if not self.values:
            return 0.0
        return sum(bool(gate(v)) for v in self.values) / len(self.values)

    def is_stable(self, gate: Callable[[float], bool]) -> bool:
        """True when every sample gives the gate the same answer.

        An unstable verdict is not a near miss -- it is a coin flip on which sample was scored,
        and reporting it as a verdict about the design is the error this module exists to stop.
        """
        return self.agreement(gate) in (0.0, 1.0)

    def verdict(self, gate: Callable[[float], bool]) -> str:
        if not self.values:
            return "unmeasured"
        return {1.0: "keep", 0.0: "drop"}.get(self.agreement(gate), "unstable")


@dataclass(frozen=True)
class StabilityReport:
    """How much of a gate's output across a set of designs is sample-dependent."""

    keep: tuple[str, ...]
    drop: tuple[str, ...]
    unstable: tuple[str, ...]

    @property
    def total(self) -> int:
        return len(self.keep) + len(self.drop) + len(self.unstable)

    @property
    def unstable_fraction(self) -> float:
        return len(self.unstable) / self.total if self.total else 0.0

    @property
    def trustworthy(self) -> bool:
        """A gate deciding more than one design in ten by coin flip is not reporting on designs.

        Deliberately strict. The campaign that produced this had 10 of 16 verdicts flip, and
        every count derived from that gate -- "six of eight failed on structure" -- was a
        statement about one sample rather than about the designs.
        """
        return self.unstable_fraction <= 0.1

    def __str__(self) -> str:
        return (f"{len(self.keep)} keep, {len(self.drop)} drop, "
                f"{len(self.unstable)} sample-dependent "
                f"({self.unstable_fraction:.0%}){'' if self.trustworthy else ' — NOT trustworthy'}")


def stability(measurements: Sequence[Sampled], gate: Callable[[float], bool]) -> StabilityReport:
    """Apply a gate to every sample of every design and report where it cannot decide.

    Run this before quoting a count from a gate. A funnel built on an unstable gate reports the
    sample that was scored, not the designs that were rejected.
    """
    keep, drop, unstable = [], [], []
    for m in measurements:
        {"keep": keep, "drop": drop}.get(m.verdict(gate), unstable).append(m.name)
    return StabilityReport(tuple(keep), tuple(drop), tuple(unstable))
