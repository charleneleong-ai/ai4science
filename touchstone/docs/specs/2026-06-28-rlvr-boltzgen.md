# RLVR for BoltzGen — fine-tuning against the touchstone reward

**Goal:** make BoltzGen propose better *metal* binders by fine-tuning it against touchstone's
verdict as the reward. touchstone's geometry prior is the CSD metal pull (`csd_reference.json`),
so **CSD metal knowledge enters through the verifier, not as raw training data** — which is the
only sound way to inject it into a *protein* generator (the CSD is small-molecule crystallography;
BoltzGen trains on the PDB).

## Why RAFT, not DPO
BoltzGen ships a PyTorch-Lightning `Training` task (`boltzgen.task.train`) with checkpoint
**resume** — but no DPO / reward / preference machinery (only the metric `filter`). So the loop is
**RAFT** (reward-ranked / reject-sampling fine-tuning): generate a pool, keep the verifier's
winners, and *supervised-fine-tune on the winners*, resuming from the released checkpoint. Repeat.
No custom RL code needed — the reward is expressed as **which samples enter the training set**.

## The loop
1. **Generate** a pool with BoltzGen (`boltzgen run`, theozyme-motif spec) → `refold_cif/` + `fold_out_npz/`.
2. **Select** the winners with touchstone (CSD-calibrated reward):
   `scripts/rlvr_select.py --npz-dir … --cif-dir … --out round_k --metal Ni2+ --keep trust`
   → `round_k/dataset/*.cif` (the kept designs) + `round_k/rewards.jsonl`. Reward =
   `reward_from_result` (consensus-weighted mean over geometry · bond-valence · nVECSUM · shape · …).
3. **Format** the kept CIFs into BoltzGen's training manifest (its data-pipeline format — the one
   step that needs BoltzGen-internal knowledge; mirror `resources/config/train` data config).
4. **Fine-tune** on the winners, resuming from the released weights, on the A100 (`bg` env):
   `python -m boltzgen.task.train.train ... --config-name boltzgen.yaml resume=<released.ckpt> data=<winners>`
5. **Repeat** — the policy drifts toward verifier-passing (CSD-plausible) metal sites.

This is exactly the verifier slot the team's `loop_runner` already exposes (its per-loop
"verifier evaluation" adapter) — touchstone is that adapter; rlvr_select is the batch version.

## Honest caveats (read before burning A100 hours)
- **The pool must contain winners.** On the current bare Ni pool *every* design is `DEFER`
  (touchstone), so strict `--keep trust` selects 0 — you can't bootstrap RAFT from it. Use the
  theozyme-motif pool (3/12 TRUST) or soften to `--keep N` (top-N by reward) for the first round,
  and improve the generation spec before tightening.
- **Reward is not yet wet-lab-calibrated.** RAFT will optimize "CSD/physics-plausible," not measured
  Kd/selectivity — so it sharpens the *prior*, not proven binding, until calibration data exists.
- **Reward-hacking / mode collapse.** Fine-tuning to a fixed verifier can collapse diversity or
  exploit verifier blind spots. Keep a diversity term, cap rounds, and hold out an
  independent check (e.g. a co-fold cross-check the reward didn't use) to detect drift.
- **Data formatting** (step 3) is the real engineering gap; the rest is wired.

## Status
- ✅ Reward + selection (`reward.py`, `scripts/rlvr_select.py`) — runnable; validated on real BoltzGen Ni output.
- ✅ BoltzGen fine-tune entrypoint identified (`Training`, resume).
- ⬜ Winners → BoltzGen training manifest (step 3).
- ⬜ First fine-tune round on the A100, then re-verify the new pool with touchstone (did TRUST-rate rise?).
