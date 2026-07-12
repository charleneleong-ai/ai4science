"""Turn the verifier stack into a single scalar reward — and rank/select designs by it.

The verifier-as-reward step toward RLVR (`biometal_rlvr`): a design's reward is the
mean score of the verifiers that ran, scaled by the consensus (trust=1, weak=0.5,
defer=0 — defense-in-depth, so a single defer collapses the reward). `rank_structures`
scores a batch best-first; `best_of_n` keeps the top design(s) — the simplest
reward-guided selection (rejection sampling) before any policy-gradient RL.
"""

from __future__ import annotations

from statistics import fmean

from .core import CONSENSUS_WEIGHT
from .service import AUTO, mlip_backbone, verify_structure


def reward_from_result(result: dict) -> float:
    """Scalar reward in [0, 1] from a verify_structure result: mean verifier score
    scaled by the consensus weight (a single defer zeroes it)."""
    scores = [v["score"] for v in result.get("verifiers", {}).values() if "score" in v]
    if not scores:
        return 0.0
    return round(CONSENSUS_WEIGHT[result["consensus"]] * fmean(scores), 4)


def rank_structures(
    structures, metal: str = "Ni2+", deep: bool = False, gate_defer: bool = False, calc=AUTO,
    precedent: bool = True, precedent_search=None, selectivity_metals=None,
) -> list[dict]:
    """Verify each structure and return results (each with a `reward`) sorted best-first.
    A structure that can't be parsed/verified scores 0 with an `error` recorded. `calc` lets a
    caller share one MLIP backbone (leave default to build per call). `precedent_search` and
    `selectivity_metals` fold the precedent and MLIP metal-swap ΔE selectivity tiers into the reward
    (see `verify_structure`) — selectivity is the lever for *target-selective* binders. With
    `gate_defer`, run the cheap geometry pass first and spend the deep MLIP tiers only on designs that
    don't already defer on geometry — a geometry-defer scores reward 0 regardless (see
    `CONSENSUS_WEIGHT`), and MLIP can only add defers, so gating can't change the selection, only save
    GPU."""
    if calc is AUTO:
        calc = mlip_backbone() if deep else None  # build the MLIP backbone once, share across the batch
    # the deep tiers are inert without `deep`, so one kwarg set serves every pass
    tiers = dict(precedent=precedent, precedent_search=precedent_search, selectivity_metals=selectivity_metals)
    scored: list[dict] = []
    for s in structures:
        try:
            if gate_defer and calc is not None:
                result = verify_structure(s, metal, **tiers)  # cheap pass first
                if result["consensus"] != "defer":
                    result = verify_structure(s, metal, deep=True, calc=calc, **tiers)
            else:
                result = verify_structure(s, metal, deep, calc=calc, **tiers)
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
