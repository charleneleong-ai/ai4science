"""Independent checks on a designed binder.

Each implements touchstone's `Verifier` protocol -- `verify(design) -> Verdict` -- so they
compose, and so the repo shares one verdict vocabulary rather than a parallel per-project one.

`Verdict.defer` means *cannot judge* (no database, unparseable structure, failed predictor),
which is deliberately distinct from judging a design and failing it. Conflating the two is how
a missing pKa silently became a passing design.

The order below is the order that matters in practice: the cheap sequence checks eliminate most
candidates, and the expensive structural ones only run on survivors.
"""

from tnf_binder.verifiers.composition import CompositionVerifier
from tnf_binder.verifiers.fold_novelty import FoldNoveltyVerifier, FoldHit
from tnf_binder.verifiers.interface import InterfaceVerifier, InterfaceMetrics
from tnf_binder.verifiers.pka_switch import PkaSwitchVerifier, HistidineSite, pka_window
from tnf_binder.verifiers.sequence_novelty import SequenceNoveltyVerifier

__all__ = [
    "CompositionVerifier",
    "SequenceNoveltyVerifier",
    "FoldNoveltyVerifier",
    "FoldHit",
    "InterfaceVerifier",
    "InterfaceMetrics",
    "PkaSwitchVerifier",
    "HistidineSite",
    "pka_window",
]
