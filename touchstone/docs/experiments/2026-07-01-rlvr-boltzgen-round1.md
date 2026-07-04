# RLVR rounds 1–4 — does the touchstone reward improve BoltzGen? (Ni, 2026-07-01)

First end-to-end run of the loop in [`docs/specs/2026-06-28-rlvr-boltzgen.md`](../specs/2026-06-28-rlvr-boltzgen.md):
generate → [`rlvr_select`](../../scripts/rlvr_select.py) → [`winners_to_targets`](../../scripts/winners_to_targets.py)
→ fine-tune → re-verify. All on `pi-a100-80gb`, Ni²⁺, the `ni_motif` theozyme spec.

## Hypothesis
Fine-tuning BoltzGen on its own touchstone-TRUST designs (RAFT) raises the fraction of a
fresh pool that touchstone trusts — i.e. the verifier-as-reward measurably improves the generator.

## Setup
- **Baseline pools:** 4 × `boltzgen run` (`ni_motif`, protein-anything, design→inverse_fold→fold),
  288 designs total, scored with `rlvr_select --keep trust` (full v0.0.3 stack: geometry z-score ·
  bond-valence · nVECSUM · polyhedron-RMSD).
- **Winners → targets:** `winners_to_targets` → 17 training targets (`structures/`+`records/`+`manifest.json`).
- **Fine-tune:** resume from released `boltzgen1_diverse.ckpt`, 300 steps, `max_lr=5e-5`,
  `lr_warmup_no_steps=30`, `use_msa=false`, `validate_structure=false`, 1×A100 bf16 (~73 min).
- **Eval:** fresh 96-design pool from the fine-tuned checkpoint (`--design_checkpoints last.ckpt`),
  re-scored with the same stack.

## Result — the reward lifts the geometry-TRUST rate 3.7×

| pool | designs | TRUST | rate |
|---|---|---|---|
| base model (R1–R4) | 288 | 17 | **5.9%** |
| **fine-tuned** | 96 | 21 | **21.9%** |

3.7× lift, ≈6.6σ under the baseline rate. Not mode-collapse: the whole distribution shifted —
`defer` fell from ~47% (base pools) to 26%, `weak` rose, mean reward 0.07→0.31, top rewards to 0.87.

## Reality check — independent MLIP physics (the reward never used)
Re-scored the 21 fine-tuned winners with the `--deep` MLIP tier (MACE relaxation + 300 K MD),
independent of the geometry reward:

| MLIP relaxation | trust | weak | collapse |
|---|---|---|---|
| 21 winners | 4 | 14 | 3 |

- **Sites are real:** 18/21 hold their donor shell under relaxation, ΔE_bind favorable for all 21
  (−3.2 to −5.5 eV) — the metal stays bound. Only 3 structurally collapse.
- **But the reward is more permissive than physics:** it called all 21 TRUST; MLIP promotes only 4
  to trust and downgrades the rest to weak (moderate drift / one lost donor). Under 300 K MD only
  `ni_motif_37` and `ni_motif_60` survive both (drift 0.31–0.35 Å, ≥98% retention).

Round 1 lifted geometry-TRUST 3.7× — real *in geometry space* — but the deep check exposed the
reward's blind spot: geometry-plausible ≠ MLIP-stable. It optimized exactly what it measured; the
fix is to measure more (rounds 2–3 below).

## Round 2 — folding MLIP into the reward (the objective tradeoff)
Added `--deep` to `rlvr_select` so the MLIP relax+MD tiers enter the consensus + `reward_from_result`
(shared MACE backbone; MLIP spent only on geometry-plausible designs). Re-selected the round-1
fine-tuned pool (96) → top-15 by the MLIP-aware reward, fine-tuned round-1's checkpoint on them
(iterative RAFT), and deep-verified a fresh pool.

The round-2 training set was genuinely MLIP-stable (12/15 MLIP-MD-trust, 0 collapse) — but mostly only
geometry-*weak* (3 trust · 12 weak), because MLIP-stability dominated the ranking. The result, on a
fresh 96-design pool:

| metric | round-1 (geometry reward) | round-2 (MLIP-aware reward) |
|---|---|---|
| geometry-TRUST | 21.9% (21/96) | **6.3% (6/96)** ↓ 3.5× |
| MLIP-MD trust | ~4/96 (subset) | **34.4% (33/96)** ↑ ~8× |

**RLVR moved exactly what the reward measured — and regressed the objective it stopped weighting.**
Dynamic stability jumped ~8×; geometry-TRUST collapsed, because the MLIP-dominated selection trained
the model on geometrically-mediocre (but stable) sites. Two rounds now demonstrate the tradeoff
empirically: a single-objective-dominated reward sacrifices the de-emphasized objective.

## Round 3 — balanced reward (the resolution)
Branched from **round-1's checkpoint** (best geometry) and fine-tuned on **8 dual-passers** —
designs that pass geometry-TRUST **and** MLIP-MD-trust, accumulated across the ft2+ft3 pools. Then
generated + deep-verified a fresh 96-design pool. The balanced set pushed *both* axes up at once:

| metric (96-design pool) | round 1 (geometry) | round 2 (MLIP) | **round 3 (balanced)** |
|---|---|---|---|
| geometry-TRUST | 21.9% | 6.3% | **66.7% (64/96)** |
| MLIP-MD trust | ~4% | 34.4% | **39% (37/96)** |
| full 6-tier TRUST (all verifiers) | ~2% | ~2% | **12.5% (12/96)** |

Geometry-TRUST hit its highest of any round (3× round-1, 11× baseline) with **no** MLIP regression,
and the "passes *every* verifier" rate (z-score · BVS · nVECSUM · polyhedron · MLIP relax · MLIP-MD)
jumped ~6×. Training on designs that were uniformly high-quality on both axes improved both — the
inverse of round-2, which trained on geometry-weak (but stable) sites and dragged geometry down.

## Verdict (rounds 1–3)
Verifier-as-reward **works, is steerable, and is reward-shaped** — RLVR optimizes precisely what the
reward selects for. The three rounds map the consequence cleanly:
- **geometry-only** → geometry ↑, MD-labile sites;
- **MLIP-only** → MLIP ↑, geometry collapses;
- **balanced (require both)** → *both* rise together, and the fully-verified rate climbs ~6×.

The lesson: **the reward must encode every objective you care about** — a single-objective-dominated
reward sacrifices whatever it stops weighting. Balanced dual-pass selection is the resolution — and
then it **saturates**: round 4 (below) holds the gains but doesn't climb further, so beyond this the
lever is a richer generation spec / new metal / wet-lab calibration, not more RAFT on the same reward.

## Round 4 — more dual-passers → plateau
Iterative RAFT continued: fine-tuned **round-3's** checkpoint on the accumulated **20 dual-passers**
(round-3's 8 + the 12 new full-6-tier winners from round-3's own pool), then generated + deep-verified
a fresh 96-design pool.

| metric (96-design pool) | round 3 (8 dual-passers) | round 4 (20 dual-passers) |
|---|---|---|
| geometry-TRUST | 66.7% | 72.9% (70/96) |
| MLIP-MD trust | 39% | 34% (33/96) |
| full 6-tier TRUST | 12.5% | 10.4% (10/96) |

**The gains plateau.** Geometry-TRUST holds high (edged up, within noise), but the fully-verified
6-tier rate is flat at ~10–12% (12→10/96 is ≈2 designs, well inside binomial noise for n=96). A second
balanced round neither climbs past round-3 nor collapses — the balanced reward's win **saturates** once
the policy already sits at the dual-pass ceiling this generation spec / motif pool can reach. Pushing
higher likely needs a *different* lever (a richer motif pool, a metal beyond Ni, or wet-lab-calibrated
reward signal), not more RAFT rounds on the same objective.

## Caveats
- **Not wet-lab-calibrated** — optimizes CSD/physics-plausibility, not measured Kd/selectivity.
- **Single seed per round**; effects are large but unreplicated.
- **Small balanced sets** (8→20 dual-passers) — dual-passers are rare (~2–12 per 96-pool), so the
  fine-tuning sets are small and overfitting is a live risk; the fresh-pool re-verify is the guard.
