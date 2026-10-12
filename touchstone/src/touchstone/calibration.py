"""Does a proposed gate actually separate designs that worked from designs that did not?

Every verifier in this package carries a threshold, and until now every one of them was a
number someone chose. Four modules even say "calibrate thresholds against wet-lab data" in
their docstrings, and none of them did, because there was nothing to calibrate against.

Proteinbase publishes an export of designs with their measured outcomes -- whether each
expressed, whether it bound, and its KD where one was fitted. That is the missing control set.
This module turns it into the question a verifier should have to answer before anyone trusts
it: *of the designs with a known outcome, does this gate keep the ones that worked?*

The three failures this was written after all look the same in hindsight, and all three would
have been caught here in one call:

  - a clash gate fixed at zero, when seven of the thirteen designs it was being compared
    against carried clashes and one carried three;
  - an expression verifier whose shipped defaults pass every design in a real set, so its
    verdict carries no information;
  - a cysteine liability argued from first principles, worth 0.7 percentage points of measured
    expression once a control set existed.

**Pool targets and the answer can invert.** Measured over 2214 designs with a Boltz-2 ipTM
prediction and a binding outcome, ipTM looks useless or backwards: designs predicted at
0.1-0.2 bound 45.7% of the time and those above 0.7 bound 33.1%. Split by target it predicts
correctly -- on Nipah glycoprotein G, 1.1% of designs below 0.5 bound against 10.6% above, a
tenfold separation. The pooled view is Simpson's paradox, because base rates run from 0.9% on
15-PGDH to 84.3% on SARS-CoV-2 RBD and the high-ipTM bins are dominated by the hardest target.
`calibrate_by_target` exists so this is the default way to look, and `calibrate` without a
target should be read as a pooled summary that may be meaningless.

A gate can also discriminate for the wrong reason, and this module will not catch that for
you. "Even or no cysteine" separates outcomes at z=4.2 on the export -- but only because it
sweeps the zero-cysteine designs into the kept group; odd against even *non-zero* is 86.7%
against 86.7%, exactly nothing. A significant separation says the two groups differ, not that
the feature you named is what differs, so the comparison still has to be the one you meant.

`discriminates` is deliberately conservative. A gate separates outcomes only when the gap
between what it keeps and what it drops survives the sampling noise in the two groups, so a
gate that looks good on eleven designs is reported as undecided rather than as working.

Provenance of every number quoted here and in the README -- the export is a snapshot, not a
stable dataset, and a later one will give different counts:

    source   proteinbase.com, full export, "Download" in the site header
    file     proteinbase_all_data_11_10_2026.csv
    rows     4069 designs, 4067 with an experimental expression verdict
    sha256   15be8539eee642b27bf66298153f9878...
    read     2026-10-12

The export is not committed: it is 32 MB and it is Proteinbase's to distribute. Nothing in the
test suite depends on it -- the tests carry the measured counts as fixtures, so they pass
without it and still fail if the arithmetic changes.
"""

from __future__ import annotations

import csv
import json
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

EXPRESSED, BINDS = "expressed", "binds"
# Two-sided 95%: a gate has to clear this to be called discriminating rather than undecided.
Z_95 = 1.96


@dataclass(frozen=True)
class Outcome:
    """One design and what the lab measured, with None wherever nothing was measured."""

    sequence: str
    expressed: bool | None = None
    binds: bool | None = None
    kd_nm: float | None = None
    target: str = ""
    design_method: str = ""

    def measured(self, outcome: str) -> bool | None:
        return self.expressed if outcome == EXPRESSED else self.binds


def _truthy(value: object) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ("true", "yes", "1"):
        return True
    if isinstance(value, str) and value.strip().lower() in ("false", "no", "0"):
        return False
    return None


def _number(value: object) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def load_export(path: str | Path, *, target: str | None = None) -> list[Outcome]:
    """Read a Proteinbase CSV export into outcomes, keeping only designs with a sequence.

    `target` filters to one target slug, which is what a new campaign wants: the base rates
    differ enough between targets that pooling them flatters or punishes a gate depending on
    which target happens to dominate the export.

    A design is counted as expressing, or binding, if any replicate says so. Experimental rows
    only -- the export also carries computational predictions under the same metric names, and
    calibrating a predictor against another predictor measures agreement, not truth.
    """
    outcomes: list[Outcome] = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sequence = (row.get("sequence") or "").strip().upper()
            if not sequence:
                continue
            try:
                evaluations = json.loads(row.get("evaluations") or "[]")
            except (TypeError, ValueError):
                continue
            seen: dict[str, list] = {}
            hit = target is None
            for evaluation in evaluations:
                if evaluation.get("type") != "experimental":
                    continue
                if target is not None and evaluation.get("target") == target:
                    hit = True
                seen.setdefault(evaluation.get("metric", ""), []).append(evaluation.get("value"))
            if not hit:
                continue
            expressed = [c for c in map(_truthy, seen.get("expressed", [])) if c is not None]
            binds = [c for c in map(_truthy, seen.get("binding", [])) if c is not None]
            kds = [k for k in map(_number, seen.get("kd", [])) if k is not None]
            outcomes.append(Outcome(
                sequence=sequence,
                expressed=any(expressed) if expressed else None,
                binds=any(binds) if binds else None,
                kd_nm=min(kds) if kds else None,
                target=target or "",
                design_method=(row.get("designMethod") or "").strip(),
            ))
    return outcomes


@dataclass(frozen=True)
class Calibration:
    """What a gate would have done to designs whose outcome is already known."""

    kept: int
    dropped: int
    kept_successes: int
    dropped_successes: int
    outcome: str = EXPRESSED

    @property
    def total(self) -> int:
        return self.kept + self.dropped

    @property
    def kept_rate(self) -> float:
        return self.kept_successes / self.kept if self.kept else 0.0

    @property
    def dropped_rate(self) -> float:
        return self.dropped_successes / self.dropped if self.dropped else 0.0

    @property
    def base_rate(self) -> float:
        return (self.kept_successes + self.dropped_successes) / self.total if self.total else 0.0

    @property
    def separation(self) -> float:
        """How much better the kept designs do. Negative means the gate keeps the wrong ones."""
        return self.kept_rate - self.dropped_rate

    @property
    def lift(self) -> float:
        """What the gate buys over not gating at all -- usually far smaller than `separation`."""
        return self.kept_rate - self.base_rate

    @property
    def standard_error(self) -> float:
        if not (self.kept and self.dropped):
            return float("inf")
        a = self.kept_rate * (1 - self.kept_rate) / self.kept
        b = self.dropped_rate * (1 - self.dropped_rate) / self.dropped
        return math.sqrt(a + b)

    @property
    def z(self) -> float:
        se = self.standard_error
        return 0.0 if se in (0.0, float("inf")) else self.separation / se

    @property
    def discriminates(self) -> bool:
        """True only when the separation survives the noise in both groups.

        A gate that keeps everything, or drops everything, discriminates nothing however good
        the kept rate looks -- which is the failure mode that let an expression verifier pass
        thirteen of thirteen designs and read as if it had judged them.
        """
        return bool(self.kept and self.dropped) and abs(self.z) >= Z_95

    @property
    def verdict(self) -> str:
        if not (self.kept and self.dropped):
            return "vacuous"
        if not self.discriminates:
            return "undecided"
        return "discriminates" if self.separation > 0 else "inverted"

    def __str__(self) -> str:
        return (f"{self.verdict}: keeps {self.kept} at {self.kept_rate:.1%}, "
                f"drops {self.dropped} at {self.dropped_rate:.1%} "
                f"(base {self.base_rate:.1%}, separation {self.separation:+.1%}, z={self.z:.1f})")


def calibrate(outcomes: Iterable[Outcome], gate: Callable[[Outcome], bool], *,
              outcome: str = EXPRESSED) -> Calibration:
    """Apply a gate to designs with a known outcome and report what it would have done.

    `gate` returns True for a design the gate keeps. Designs with nothing measured for
    `outcome` are excluded rather than counted as failures: an unmeasured design is not a
    failed one, and collapsing the two is how a gate acquires a track record it never earned.
    """
    kept = dropped = kept_ok = dropped_ok = 0
    for design in outcomes:
        measured = design.measured(outcome)
        if measured is None:
            continue
        if gate(design):
            kept += 1
            kept_ok += measured
        else:
            dropped += 1
            dropped_ok += measured
    return Calibration(kept, dropped, kept_ok, dropped_ok, outcome)


@dataclass(frozen=True)
class Band:
    """One slice of a continuous feature, and how designs in it did."""

    low: float
    high: float
    n: int
    successes: int

    @property
    def rate(self) -> float:
        return self.successes / self.n if self.n else 0.0

    def __str__(self) -> str:
        return f"{self.low:g} to {self.high:g}: {self.rate:.1%} (n={self.n})"


def bands(outcomes: Iterable[Outcome], feature: Callable[[Outcome], float],
          edges: Sequence[float], *, outcome: str = EXPRESSED) -> list[Band]:
    """Outcome rate per slice of a continuous feature, before choosing any threshold.

    Worth looking at first: a feature whose rate is flat across its whole range has no
    threshold worth choosing, and a monotonic one tells you where to put the gate rather than
    leaving it to be guessed.
    """
    measured = [(feature(d), d.measured(outcome)) for d in outcomes]
    measured = [(v, bool(m)) for v, m in measured if m is not None]
    out = []
    for low, high in zip(edges, edges[1:]):
        inside = [m for v, m in measured if low <= v < high]
        out.append(Band(low, high, len(inside), sum(inside)))
    return out


def calibrate_by_target(outcomes: Iterable[Outcome], gate: Callable[[Outcome], bool], *,
                        outcome: str = EXPRESSED, min_group: int = 15,
                        ) -> dict[str, Calibration]:
    """One calibration per target, which is the only safe way to read a gate across campaigns.

    Base rates differ by two orders of magnitude between targets, so a gate pooled over all of
    them can invert its own sign. Targets where either side of the gate holds fewer than
    `min_group` designs are omitted rather than reported at a rate that cannot be distinguished
    from noise.
    """
    grouped: dict[str, list[Outcome]] = {}
    for design in outcomes:
        if design.measured(outcome) is not None:
            grouped.setdefault(design.target, []).append(design)
    per_target = {}
    for target, designs in grouped.items():
        result = calibrate(designs, gate, outcome=outcome)
        if min(result.kept, result.dropped) >= min_group:
            per_target[target] = result
    return per_target
