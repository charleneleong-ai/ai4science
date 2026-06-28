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
3. **Format** the kept CIFs into BoltzGen training targets with `scripts/winners_to_targets.py` (validated — see below).
4. **Fine-tune** on the winners, resuming from the released weights, on the A100 (`bg` env):
   `python -m boltzgen.task.train.train ... --config-name boltzgen.yaml resume=<released.ckpt> data=<winners>`
5. **Repeat** — the policy drifts toward verifier-passing (CSD-plausible) metal sites.

This is exactly the verifier slot the team's `loop_runner` already exposes (its per-loop
"verifier evaluation" adapter) — touchstone is that adapter; rlvr_select is the batch version.

## Winners → BoltzGen training targets (step 3 — `scripts/winners_to_targets.py`, validated)
BoltzGen's trainer (`task/train/data.py`) reads, per design, a **`Record` JSON + a `Structure`
`.npz`** under `DatasetConfig.target_dir`, plus a `Manifest` at `manifest_path`. The validated
layout the converter writes (confirmed by a dump→load round-trip on real winner CIFs):

    out/structures/<id>.npz   # Structure.dump  (parse_mmcif → ParsedStructure.data)
    out/records/<id>.json     # Record.dump     (StructureInfo + per-chain ChainInfo)
    out/manifest.json         # Manifest(records=[…])  → DatasetConfig.manifest_path

Path, in the `bg` env (`boltzgen.data`): `parse_mmcif(cif, moldir=~/.boltz/mols) -> ParsedStructure`
→ `.data` is the `Structure` (dump it), `.info` is the `StructureInfo` (goes in the `Record`). Build
one `ChainInfo` per `structure.chains` row (protein `mol_type=0` + metal ligand `mol_type=3`).

**No MSA needed.** De-novo designs carry no alignment; the trainer gates MSA loading on `use_msa`,
so fine-tune with **`use_msa: false`** (single-sequence) and `msa_dir` is never read — `msa_id=-1`.
The generation `.npz` is **not** reusable as a training `Structure` (no `chains` array — verified),
so the converter re-parses the winner CIFs. Needs the `~/.boltz/mols` CCD cache (`NI.pkl` etc.) to
resolve components.

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
- ✅ Step-3 converter `scripts/winners_to_targets.py` — built + validated in the `bg` env (dump→load round-trip on real round-1 winners; `structures/`+`records/`+`manifest.json`).
- ✅ Round-1 motif pool (48) generated → `rlvr_select` kept **4 TRUST** winners (`ni_motif_13/16/40/43`).
- 🟡 Round-2 motif pool (48) generating to enlarge the winner set before fine-tuning.
- ⬜ Merge round-1+round-2 winners → run the converter → first resume-train round on the A100, then re-verify the new pool with touchstone (did TRUST-rate rise?).
