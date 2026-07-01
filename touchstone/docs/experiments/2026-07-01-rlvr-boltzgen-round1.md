# RLVR rounds 1–2 — does the touchstone reward improve BoltzGen? (Ni, 2026-07-01)

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

## Verdict
RLVR **works** — the generator genuinely shifted toward the reward's notion of a good site, and the
3.7× lift is real *in geometry space*. But geometry-plausible ≠ MLIP-stable: the deep check exposes
the reward's blind spot (it doesn't penalize dynamic lability). RAFT optimized exactly what we
asked; the fix is to ask for more.

## Caveats
- **Not wet-lab-calibrated** — optimizes CSD/physics-plausibility, not measured Kd/selectivity.
- **No matched base-model deep baseline** — can't yet claim the fine-tune made MLIP-stability
  *worse*, only that geometry-TRUST overstates it. A base-model `--deep` pass would close that.
- **Single run / single seed.** Effect is large but unreplicated.

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

## Verdict (both rounds)
Verifier-as-reward **works and is steerable** — it reliably optimizes whatever the reward encodes.
That is also its hazard: geometry-only → MD-labile designs; MLIP-only → geometry-poor designs. The
reward must be **balanced**, and selection must require *both* (geometry-TRUST ∧ MLIP-stable), not let
either term dominate.

## Next move — round 3, balanced reward
Branch from **round-1's checkpoint** (best geometry) and train only on **dual-passers** (geometry-TRUST
∧ MLIP-MD-trust). These are rare (~2–4 per 96-pool), so the first balanced set is small (8 accumulated
across the ft2+ft3 pools). Success = a pool that holds geometry-TRUST **and** MLIP-stability together —
neither collapsing. If the dual-pass rate is too low to train well, grow the set with more generation
or adopt a weighted continuous reward that keeps both terms.
