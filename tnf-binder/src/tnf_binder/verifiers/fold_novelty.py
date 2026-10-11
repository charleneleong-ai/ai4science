"""Score fold novelty the way Proteinbase does, because fold is the binding constraint.

Two submitted designs differed by 0.6% in sequence identity (37.5% vs 38.1%). One passed
Proteinbase's automatic novelty check and one failed. The failure was a beta-sandwich matching
CATH 2.60.40 -- the immunoglobulin-like superfamily -- at TM-score 0.893 over 89% coverage.
Sequence novelty could not see it; fold similarity decided it.

Proteinbase's published scale, which submissions must reach 3 on:

  4  de novo   sequence identity <=30%  AND  less than moderate structural similarity
  3            identity >30%  XOR  moderate structural similarity
  2            identity >70%, or high structural similarity,
               or identity >30% together with moderate structural similarity
  1  known     identity >70%  AND  at least moderate structural similarity

Structural similarity is TM-score weighted by coverage:
  high      >=70% coverage at TM >= 0.8
  moderate   >70% coverage at TM >= 0.5

CALIBRATION, honestly stated. Checked against three designs with known Proteinbase verdicts:

  the FAILING design is predicted correctly -- high fold similarity (TM 0.893 over 89%) gives
  score 2, and it failed.

  the two PASSING designs are UNDER-scored. At ~38% identity with locally-measured moderate
  similarity this returns 2, yet Proteinbase passed them. Their pipeline is ESMFold ->
  FoldSeek -> TM-align with domain segmentation over PDB/patent/antibody databases; this runs
  FoldSeek alntmscore over CATH50 on co-folded structures, and measures more similarity than
  they do.

So this is a CONSERVATIVE pre-filter: a `high` verdict reliably predicts rejection, while a
`moderate` verdict may still pass upstream. Treat a failure here as decisive and a borderline
score as a warning, not a disqualification. An earlier claim that it reproduced all three
verdicts was wrong -- the two passes were scored correctly only because a missing identity
value defaulted to zero.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from touchstone.core import Verdict
from touchstone.search import FoldHit

from tnf_binder.design import TnfDesign

HIGH_TM, HIGH_COVERAGE = 0.8, 0.70
MODERATE_TM, MODERATE_COVERAGE = 0.5, 0.70
IDENTITY_FAMILIAR, IDENTITY_KNOWN = 30.0, 70.0


def structural_similarity(hits: list[FoldHit]) -> tuple[str, FoldHit | None]:
    """Proteinbase's coverage-weighted verdict: none | moderate | high."""
    high = [h for h in hits if h.tm_score >= HIGH_TM and h.coverage >= HIGH_COVERAGE]
    if high:
        return "high", max(high, key=lambda h: h.tm_score)
    moderate = [h for h in hits if h.tm_score >= MODERATE_TM and h.coverage > MODERATE_COVERAGE]
    if moderate:
        return "moderate", max(moderate, key=lambda h: h.tm_score)
    return "none", max(hits, key=lambda h: h.tm_score) if hits else None


def novelty_score(identity: float, similarity: str) -> int:
    """The published 1-4 scale. 4 is de novo; submissions need >=3."""
    is_high = similarity == "high"
    is_moderate = similarity in ("moderate", "high")
    if identity > IDENTITY_KNOWN and is_moderate:
        return 1
    if identity > IDENTITY_KNOWN or is_high or (identity > IDENTITY_FAMILIAR and is_moderate):
        return 2
    if identity > IDENTITY_FAMILIAR or is_moderate:
        return 3
    return 4


class FoldNoveltyVerifier:
    """Trusts a design whose predicted novelty score reaches `min_score`.

    `fold_search` returns structural hits for a design; `identity_of` returns its best
    sequence identity as a percentage. Both are injected so the databases and their snapshot
    dates stay an explicit choice.
    """

    def __init__(
        self,
        fold_search: Callable[[TnfDesign], list[FoldHit]],
        identity_of: Callable[[TnfDesign], float],
        *,
        min_score: int = 3,
    ) -> None:
        self.fold_search = fold_search
        self.identity_of = identity_of
        self.min_score = min_score

    def verify(self, design: TnfDesign) -> Verdict:
        if design.structure is None:
            return Verdict.defer("no structure to search with")
        try:
            hits = self.fold_search(design)
            identity = self.identity_of(design)
        except Exception as exc:
            return Verdict.defer(f"novelty search failed: {type(exc).__name__}")

        similarity, closest = structural_similarity(hits)
        score = novelty_score(identity, similarity)
        metrics = {
            "identity_pct": round(identity, 1),
            "structural": similarity,
            "novelty_score": score,
            "best_tm": round(closest.tm_score, 3) if closest else None,
            "tm_coverage": round(closest.coverage, 3) if closest else None,
            "closest_fold": closest.target if closest else None,
        }
        reason = (f"novelty {score}/4 — {identity:.0f}% identity, {similarity} fold similarity"
                  + (f" (TM {closest.tm_score:.2f} over {closest.coverage:.0%} "
                     f"to {closest.target})" if closest else ""))
        return Verdict(score / 4.0, trust=score >= self.min_score, ood=False,
                       metrics=metrics, reason=reason)
