"""Reading Boltz confidence output with the binder/target split made explicit.

**Global ipTM is unusable the moment the target is more than one chain.** Every target-internal
chain pair scores on how well the target folds against itself -- which is high, and has nothing
to do with the binder -- and those pairs dominate the global average. On a homotrimer target,
three of the six pairs in a binder/target complex are internal and sit at ~0.85-0.95 whatever
the binder does. Measured across 57 designs, global ipTM overstated the designed interface by
**+0.405**; truncating the same target from trimer to dimer once lifted global ipTM from 0.20
to 0.52 while the binder interface itself did not move at all.

So every number here comes from `pair_chains_iptm` between the binder chain and each target
chain, symmetrised (Boltz's matrix is asymmetric by up to 0.08 and neither direction is
privileged), and `inflation` is reported alongside so the artefact stays visible rather than
merely avoided.

**Aggregation across diffusion samples is by MEAN, not best-of-n.** Boltz writes one confidence
file per sample; ranking on the best sample selects for sampling variance rather than design
quality. Best and mean differ by up to 0.35 within a single design, and the ordering changes --
a design ranked 1st on best-of-3 fell mid-pack on the mean. This mirrors `boltzgen.py`, which
aggregates its own refold samples the same way and for the same reason.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

BINDER_CHAIN = "0"


def mean_or_zero(values: Iterable[float]) -> float:
    """0.0 for an empty pool, so a design with no samples ranks last instead of raising."""
    pool = tuple(values)
    return mean(pool) if pool else 0.0


@dataclass(frozen=True)
class BoltzConfidence:
    """One Boltz-2 confidence JSON (a single diffusion sample)."""

    name: str
    global_iptm: float
    chain_ptm: dict[str, float]
    pair_iptm: dict[str, dict[str, float]]
    binder_chain: str = BINDER_CHAIN
    # Predicted distance error across the interface, lower being better. Boltz does not always
    # write it, so it defaults to NaN rather than to a number that would read as a measurement.
    complex_ipde: float = float("nan")

    @classmethod
    def from_file(
        cls, path: Path, binder_chain: str = BINDER_CHAIN
    ) -> "BoltzConfidence":
        payload = json.loads(path.read_text())
        return cls(
            name=path.stem,
            global_iptm=payload["iptm"],
            chain_ptm=payload["chains_ptm"],
            pair_iptm=payload["pair_chains_iptm"],
            binder_chain=binder_chain,
            complex_ipde=payload.get("complex_ipde", float("nan")),
        )

    @property
    def target_chains(self) -> tuple[str, ...]:
        return tuple(c for c in self.chain_ptm if c != self.binder_chain)

    def pair(self, a: str, b: str) -> float:
        """Symmetrised pairwise ipTM.

        Boltz's matrix is asymmetric -- chain 0 vs 1 and 1 vs 0 differ by up to 0.08 -- and
        neither direction is privileged, so both are averaged.
        """
        return mean((self.pair_iptm[a][b], self.pair_iptm[b][a]))

    @property
    def binder_to_target(self) -> float:
        """The number that actually measures the designed interface."""
        return mean_or_zero(self.pair(self.binder_chain, t) for t in self.target_chains)

    @property
    def target_internal(self) -> float:
        """Confidence in the target's own trimer interfaces -- what inflates global ipTM."""
        targets = self.target_chains
        return mean_or_zero(
            self.pair(a, b) for i, a in enumerate(targets) for b in targets[i + 1 :]
        )

    @property
    def binder_ptm(self) -> float:
        """Monomer confidence in the binder's own fold, independent of docking."""
        return self.chain_ptm[self.binder_chain]


@dataclass(frozen=True)
class DesignConfidence:
    """Every diffusion sample for one design, aggregated by mean."""

    name: str
    samples: tuple[BoltzConfidence, ...]

    @classmethod
    def from_dir(
        cls, directory: Path, binder_chain: str = BINDER_CHAIN
    ) -> "DesignConfidence":
        samples = tuple(
            BoltzConfidence.from_file(p, binder_chain)
            for p in sorted(directory.glob("confidence_*.json"))
        )
        return cls(directory.name, samples)

    def __bool__(self) -> bool:
        return bool(self.samples)

    @property
    def binder_to_target(self) -> float:
        """Mean across samples -- the ranking statistic."""
        return mean_or_zero(s.binder_to_target for s in self.samples)

    @property
    def binder_to_target_best(self) -> float:
        """Best sample. Reported for comparison only; ranking on it selects for variance."""
        return max((s.binder_to_target for s in self.samples), default=0.0)

    @property
    def global_iptm(self) -> float:
        return mean_or_zero(s.global_iptm for s in self.samples)

    @property
    def target_internal(self) -> float:
        return mean_or_zero(s.target_internal for s in self.samples)

    @property
    def binder_ptm(self) -> float:
        return mean_or_zero(s.binder_ptm for s in self.samples)

    @property
    def complex_ipde(self) -> float:
        """Best (lowest) interface PDE across samples -- context, never the ranking statistic."""
        return min((s.complex_ipde for s in self.samples), default=float("nan"))

    @property
    def inflation(self) -> float:
        """How much global ipTM overstates the interface. Large means the artefact is biting."""
        return self.global_iptm - self.binder_to_target


def rank_designs(
    predictions: Path, binder_chain: str = BINDER_CHAIN
) -> list[DesignConfidence]:
    """Every folded design under a Boltz predictions directory, best mean interface first.

    A missing directory yields nothing rather than raising, so callers report it the same way
    they report a directory that holds no confidence files.
    """
    if not predictions.is_dir():
        return []
    found = (
        DesignConfidence.from_dir(d, binder_chain)
        for d in sorted(predictions.iterdir())
        if d.is_dir()
    )
    return sorted(
        (d for d in found if d), key=lambda d: d.binder_to_target, reverse=True
    )
