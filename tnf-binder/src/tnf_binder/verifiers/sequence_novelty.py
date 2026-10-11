"""Catch memorised natural sequences before they reach a submission.

Generative models return known proteins. In this campaign the single highest-scoring
construct -- best interface score, lowest interface PDE, an ideal-pKa histidine -- was
verbatim human ubiquitin, and would have been disqualified on sight.

An exact k-mer match is decisive rather than statistical: a genuine de novo sequence
essentially never shares a 20-mer with a natural protein. It is a memorisation screen, not a
homology search, and it will not see a remote structural relative -- that is what
FoldNoveltyVerifier is for.
"""

from __future__ import annotations

from collections.abc import Callable

from touchstone.core import Verdict

from tnf_binder.design import TnfDesign


class SequenceNoveltyVerifier:
    """Trusts a design that shares no exact k-mer with the reference corpus.

    `lookup` returns the k-mers of the query that occur in the corpus. It is injected so the
    corpus (and its snapshot date) is an explicit choice of the caller, not a hidden global.
    """

    def __init__(
        self,
        lookup: Callable[[list[str]], set[str]],
        *,
        k: int = 20,
        stride: int = 10,
        corpus: str = "unspecified corpus",
    ) -> None:
        self.lookup = lookup
        self.k = k
        self.stride = stride
        self.corpus = corpus

    def verify(self, design: TnfDesign) -> Verdict:
        kmers = design.kmers(self.k, self.stride)
        if not kmers:
            return Verdict.defer(f"sequence shorter than k={self.k}")
        try:
            hits = self.lookup(kmers)
        except Exception as exc:
            return Verdict.defer(f"corpus lookup failed: {type(exc).__name__}")

        matched = [km for km in kmers if km in hits]
        metrics = {"k": self.k, "kmers": len(kmers), "hits": len(matched),
                   "corpus": self.corpus}
        if matched:
            return Verdict(0.0, trust=False, ood=False, metrics=metrics,
                           reason=f"memorised: {len(matched)} exact {self.k}-mer "
                                  f"matches in {self.corpus}, e.g. {matched[0]}")
        return Verdict(1.0, trust=True, ood=False, metrics=metrics,
                       reason=f"no exact {self.k}-mer match in {self.corpus}")
