"""Turn the verifier stack into a single scalar reward — and rank/select designs by it.

The verifier-as-reward step toward RLVR (`biometal_rlvr`): a design's reward is the
mean score of the verifiers that ran, scaled by the consensus (trust=1, weak=0.5,
defer=0 — defense-in-depth, so a single defer collapses the reward). `rank_structures`
scores a batch best-first; `best_of_n` keeps the top design(s) — the simplest
reward-guided selection (rejection sampling) before any policy-gradient RL.
"""

from __future__ import annotations

from statistics import fmean

from .service import verify_structure

_CONSENSUS_WEIGHT = {"trust": 1.0, "weak": 0.5, "defer": 0.0}


def reward_from_result(result: dict) -> float:
    """Scalar reward in [0, 1] from a verify_structure result: mean verifier score
    scaled by the consensus weight (a single defer zeroes it)."""
    scores = [v["score"] for v in result.get("verifiers", {}).values() if "score" in v]
    if not scores:
        return 0.0
    return round(_CONSENSUS_WEIGHT[result["consensus"]] * fmean(scores), 4)


def rank_structures(structures, metal: str = "Ni2+", deep: bool = False) -> list[dict]:
    """Verify each structure and return results (each with a `reward`) sorted best-first.
    A structure that can't be parsed/verified scores 0 with an `error` recorded."""
    scored: list[dict] = []
    for s in structures:
        try:
            result = verify_structure(s, metal, deep)
        except Exception as e:  # unparseable / no metal ⇒ worst reward, recorded
            result = {"structure": str(s), "consensus": "defer", "verifiers": {}, "error": f"{type(e).__name__}: {e}"}
        result["reward"] = reward_from_result(result)
        scored.append(result)
    scored.sort(key=lambda r: r["reward"], reverse=True)
    return scored


def best_of_n(structures, metal: str = "Ni2+", deep: bool = False, n: int = 1) -> list[dict]:
    """Reward-guided selection: the top-n structures by verifier reward — rejection
    sampling, the tractable first step before policy-gradient RLVR."""
    return rank_structures(structures, metal, deep)[:n]
