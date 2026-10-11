"""One designed binder, and the provenance that must travel with it."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class TnfDesign:
    """A candidate binder.

    `generator_confidence` is recorded and never trusted: AlphaFold2 scored its own optimised
    designs at i_pTM 0.92 where an independent model scored the same structures 0.21.
    """

    name: str
    sequence: str
    generator: str = "unknown"
    generator_confidence: float = 0.0
    structure: Path | None = None
    his_positions: tuple[int, ...] = ()
    metrics: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.sequence:
            raise ValueError(f"{self.name}: empty sequence")

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def composition(self) -> Counter:
        return Counter(self.sequence)

    @property
    def has_nonstandard_residues(self) -> bool:
        return bool(set(self.sequence) - STANDARD_AA)

    def kmers(self, k: int, stride: int = 1) -> list[str]:
        """Overlapping k-mers, always including the final window."""
        if self.length < k:
            return []
        out = [self.sequence[i:i + k] for i in range(0, self.length - k + 1, stride)]
        tail = self.sequence[-k:]
        if tail not in out:
            out.append(tail)
        return out
