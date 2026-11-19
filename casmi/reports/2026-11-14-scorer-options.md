# CASMI 2026 scorer options (2026-11-14)

**Recommendation:** build a calibrated, structure-disjoint holdout first, then test one scorer from a different model
family, FLARE fine-tuned on `train.parquet`, as a 50/50 per-molecule z-blend on ranks 2–25 of the fork487 pipeline.
Submit only if the holdout gain clears 2 × its standard error. MIST `msg_aug` is the cheap fallback. Stop adding
near-clones of the current scorer.

This is a scorer-focused delta on the [10-30 sweep](2026-10-30-casmi-sota-delta.md). Raw notes:
[`methods.md`](../research_notes/2026-11-14-scorer-options/methods.md) (papers, models, licences) and
[`kaggle_discussion.md`](../research_notes/2026-11-14-scorer-options/kaggle_discussion.md) (forum and notebooks).
Paper numbers are the authors' own unless marked *audit*. Forum numbers are single public-LB submissions and are
not verified here.

## Where we are

- Public 0.487 (fork487 v1, submission 58272211); leaders at 0.540, about 7 public molecules ahead.
- The pipeline's neural term is a 0.5/0.5 per-molecule z-blend of two models of one architecture:
  - sanpier's `casmi-spectrum-model-v1`;
  - popwenman's 80-epoch retrain of it (`casmi26-scorer-train-v80`).

  Both are peak-set transformers (width 192, 3 layers, 6 heads), trained contrastively against 2048-bit Morgan
  fingerprints with 16 hard negatives.
- **Recall is not the bottleneck; ranking is.** The pipeline's global pool holds 99.5% of Enveda-180 train answers
  and 93.6% of `enveda-np-examples` answers (2026-11-14 LOTUS reach check). The settings layers on top (mass prior,
  documentation prior, reranker weights) all landed at 0.479–0.487.
- **The kit's own holdout does not predict the leaderboard:**
  - v80 reaches a holdout diagnostic MRR@25 of 0.660, the shipped kit 0.484;
  - yet longer training scored lower in public (40 × 2000 steps 0.453 against 0.472; a self-trained scorer replacing
    the public one 0.460);
  - our seed blend lost 0.016.

  The holdout is a random split of train identities, so near-duplicate structures leak across it.

## What the forum shows for scorer swaps (since 10-30)

From the results tables in the zz3t3rn41 / guhongbin lineage notebooks
([youyouwou v225](https://www.kaggle.com/code/youyouwou/casmi26-v225v589np)):

| Change | Public LB |
|---|---|
| 50/50 z-blend, public scorer + a self-trained retrain (bgunshui), or + popwenman v80 | 0.479 → 0.480 |
| Self-trained scorer **replacing** the public one | 0.460 (−0.019) |
| Second-scorer weight 0.3 / 0.4 / 0.6 instead of 0.5 | 0.474 each (identical; possibly a fallback) |
| Natural-product fine-tuned kit scorer as a gated rank-1 challenger ([v211](https://www.kaggle.com/code/youyouwou/casmi26-v211v589ftchallenge)) | 0.487 → 0.424 (not a clean single change) |
| Seven retrained FPNet swaps, 26.5% better offline on fingerprint Tanimoto | all lost, 0.418–0.443 |
| JESTR merged by reciprocal-rank fusion on a 0.440 base (10-17 sweep) | −0.022 |

So far, blending in a near-clone moves the score by about ±0.003, and giving any new model more weight costs 0.01–0.06.
**No one has published a leaderboard number for a scorer from a different family blended into the 0.487 base.** That
gap is the opening, and also the risk.

## Options

| # | Option | What it adds | Retrieval, MassSpecGym | Weights · licence | Rules | Cost | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | **FLARE** (Hassoun lab) | Peak-to-atom contrastive matching: scores which fragments explain which substructures, not a global fingerprint | mass track hit@1/5/20 43.2 / 75.6 / 92.9; formula track 22.7 / 50.0 / 75.2 (*audit* agrees) | [HF Space HassounLab/FLARE](https://huggingface.co/spaces/HassounLab/FLARE) · MIT · trained on MassSpecGym | MassSpecGym-trained models cleared 11-11 ([757307](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/757307)); no FLARE-specific ruling | Fine-tune 1–2 days on the A100; minutes on a T4 | **Test first** |
| 2 | **MIST `msg_aug`** | Spectrum → fingerprint, cosine to each candidate | No ranking number; our probe: fingerprint Tanimoto 0.385 on 50 timsTOF Enveda-180 spectra (old checkpoint 0.26) | [Zenodo 23303289](https://zenodo.org/records/23303289) · MIT | Allowed by name ([751373](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/751373)) | Half a day; seconds per molecule (batch size 1, attention-mask bug) | **Cheap fallback** |
| 3 | **DreaMS** + retrieval head | Self-supervised encoder pretrained on mined MassIVE/GNPS spectra; we train the molecule tower | Frozen alignment caps at MRR 0.227 on decontaminated MassSpecGym ([2610.37979](https://arxiv.org/abs/2610.37979)); no fine-tuned number | [HF roman-bushuiev/DreaMS](https://huggingface.co/roman-bushuiev/DreaMS) · MIT | Official checkpoints allowed 11-05 ([754610](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754610)); fine-tuning unanswered | 3–5 days of training | Only if 1 and 2 fail and a week remains |
| 4 | **FIORA-OS** | Forward spectrum prediction, then cosine to the query | No identification numbers for OS | [BAMeScience/fiora](https://github.com/BAMeScience/fiora) · MIT | Allowed by name (751373) | ~15 min for the same-formula top 25 | Skip: needs continuous eV, but the test gives nominal values ([753419](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753419)); −0.014 on an offline panel in negative mode |
| 5 | **MS2DeepScore 2.0** (v5) | Spectrum-to-spectrum similarity, not spectrum-to-structure | No retrieval number | [Zenodo 13897744](https://zenodo.org/records/13897744/latest) · CC BY 4.0 | Allowed by name | Seconds | Skip as a scorer; only useful for library matching across instruments |
| 6 | **JESTR** | Global contrastive scorer | formula hit@1/5/20 11.8 / 33.5 / 61.5 (*audit*) | [HassounLab/JESTR1](https://github.com/HassounLab/JESTR1) · MIT | Allowed by inference | Minutes | Skip: weaker than FLARE and already lost on the board |

**Excluded:**
- MVP and CSU-MS2: no licence file.
- MassFormer: NIST-trained.
- FraGNNet and the Separability recipe: no weights. The Separability recipe has the strongest decontaminated number
  found, MRR 0.566.
- GLMR, SpecBridge, MS-BART and DiffMS: shortcut-inflated per the audit, or de novo.
- SIRIUS / CSI:FingerID: banned.

**Why FLARE despite the JESTR loss.** The 10-17 sweep ranked FLARE low because JESTR, from the same lab, lost 0.022.
Two things differ:
- FLARE's audited formula-track hit@1 is twice JESTR's (22.7 against 11.8).
- JESTR went in by rank fusion at full weight. We would blend FLARE at 0.5 into ranks 2–25 only, leaving rank 1 to
  the pipeline, which is the form that has not lost on the board.

FLARE's value here is that it is a different family. Its errors should correlate less with the kit's than a retrain's
do, and complementarity, not standalone strength, is what a blend needs.

## Plan

**Step 0: a holdout that agrees with the leaderboard (1–2 days).**
- Hold out about 1,000 Enveda-180 molecules from train, disjoint from the training identities by InChIKey14 *and*
  with no training structure at Tanimoto ≥ 0.7. That mimics novel natural products.
- Retrain the kit on train minus the holdout with [`run_seed.py`](../scripts/seed_retrain/run_seed.py) (~1.5 h on
  the A100).
- Score pool-style candidate lists: the same formula/mass window, 25–200 candidates.
- **Calibration gate, fixed now:** before trusting the holdout for anything, it must reproduce the sign of two known
  leaderboard deltas:
  - 16 × 500 steps beats 40 × 2000 steps;
  - the public scorer alone beats the seed blend.

  If it fails, no offline number can steer a submission, and we stop scorer work.

**Step 1: FLARE (2–3 days).**
- Fine-tune FLARE from its MassSpecGym checkpoint on the same train minus holdout.
- Measure on the holdout: FLARE alone, the kit alone, and the 50/50 z-blend on ranks 2–25.
- **Submission gate, fixed now:** the blend's paired MRR@25 gain over the kit is at least 2 × SE, about +0.013 on
  1,000 molecules (10-30 sweep: SE ≈ 0.0065 on 600).
- Then:
  1. upload the weights as a private dataset;
  2. add one scoring cell to fork487 with a fallback to the current score on any failure;
  3. gate the commit run on the FLARE cell's own log line (molecules scored, nulls counted);
  4. submit one change.

**Step 2: MIST `msg_aug` (half a day).** This is the same holdout test with fingerprint cosine as the second z-term.
It runs whether or not FLARE passes, since it is cheap. If both pass, test a three-way blend as a separate submission.

**Step 3: DreaMS.** Only if Steps 1–2 fail and at least 7 days remain before 14 Dec.

**Budget:**
- 30 days and about 5 submissions a day remain, so submissions are not the constraint; holdout calibration and
  A100 time are.
- Final picks stay fork487 v1 (pick A) plus a diverse line (pick B), unless a scorer blend clears both the holdout
  gate and +0.005 on the board.

## Step 0 result (2026-11-14, gate FAIL)

The holdout ran on the A100 ([`run_calibration.sh`](../scripts/holdout/run_calibration.sh), commit
[`8cd63df`](https://github.com/charleneleong-ai/ai4science/commit/8cd63df)):
- training: 271,045 molecules;
- excluded: 1,038 within Tanimoto 0.7 of a held-out molecule;
- validation: 1,000;
- eval: 750 Enveda-180 molecules, of which 727 have same-formula pool lists.

| Scorer vs kit (16 × 500) | Holdout MRR@25 | Δ [95% CI] | Public LB |
|---|---|---|---|
| kit | 0.664 | — | 0.472 |
| long (40 × 2000) | 0.822 | **+0.158** [+0.131, +0.184] | 0.453 (−0.019) |
| kit + s1 + s2 + s3 z-blend | 0.701 | **+0.037** [+0.022, +0.053] | −0.016 (our fork479 v2) |
| s1 / s2 / s3 alone | 0.655 / 0.678 / 0.651 | −0.009 / +0.015 / −0.012 (all CIs span 0) | — |

**Both checks fail, by wide margins.** On structure-disjoint Enveda-180 molecules, longer training and seed
blending are large, significant gains; on the leaderboard both lost. Under the gate fixed beforehand, this holdout
cannot steer a submission, so FLARE and MIST are not run on its evidence.

What the disagreement does and does not say:
- The two leaderboard deltas are single submissions, each about one noise SD (±0.016). So the targets themselves
  are weak. But the holdout's gains are 10× larger than anything the board showed, so noise alone does not reconcile
  them.
- Most likely, the scorer's standalone ranking is not what the pipeline's score depends on. The neural term enters
  only on ranks 2–25 behind library matches and leader rules. And/or the hidden test is not Enveda-180-like: the
  host describes it as natural products.
- **Bug:** the natural-product group never reached eval. `build_holdout.py` labelled each molecule by its first
  `ingest_lib`. The 250 `enveda-np-examples` molecules are listed first under other libraries (145 `drug_plus`,
  88 `gnps`, 17 others), and none of them appears under `enveda-180`. So none was labelled natural product, and all
  250 stayed in training.

### Natural-product rerun (2026-11-19, post-hoc, exploratory)

The fix labels a molecule `enveda-np-examples` when any of its spectra comes from that library. The eval set becomes
the originally intended one: the same 750 Enveda-180 molecules plus all 250 natural-product molecules, with their
Tanimoto ≥ 0.7 neighbours removed from training. The five variants and the scorer stay unchanged.

The gate is the same two checks, judged on the natural-product eval molecules only (`eval_holdout.py gate
enveda-np-examples`): kit beats long, and kit beats kit + s1 + s2 + s3. It was fixed before this run, but this
subset was chosen after the first gate failed, so a pass is exploratory evidence, not a pre-registered result.

- **On a pass:** FLARE is tried on this holdout before any submission.
- **On a fail:** scorer work is closed.

**Result (2026-11-19, gate FAIL).** Eval came to 984 molecules: 750 Enveda-180 and 234 natural products, 16 of
which had no same-formula decoys. Training dropped 2,742 near neighbours. On the 226 scored natural-product molecules:

| Scorer vs kit | Holdout MRR@25 | Δ [95% CI] | Public LB |
|---|---|---|---|
| kit | 0.311 | — | 0.472 |
| long | 0.411 | **+0.100** [+0.059, +0.142] | −0.019 |
| kit + s1 + s2 + s3 | 0.328 | +0.017 [−0.007, +0.041] | −0.016 |

Across all 953 scored eval molecules, long gains +0.138 and the blend +0.030. Natural products are much harder for
the scorer alone (0.31 against 0.66 on Enveda-180), but the direction is the same as on 14 November and opposite to
the board. **A same-formula holdout does not predict the leaderboard for this scorer, on either library, so scorer
work is closed on holdout evidence.** Any further scorer family (MIST, DreaMS) can only be judged by one gated
submission against the current best.

## Not established

- Whether any offline holdout tracks this leaderboard; Step 0 tests it.
- FLARE's licence and data provenance come from the 10-17 sweep and were not re-read. No host ruling names FLARE.
  Whether its mass-track number used true or predicted formulas is unknown.
- How FLARE fine-tunes on timsTOF spectra, and whether it fits a T4 within the notebook's time limit. All T4 runtimes
  here are estimates.
- MIST `msg_aug` as a ranker: only fingerprint Tanimoto has been measured.
- Whether the host allows fine-tuning DreaMS `ssl_model.ckpt`, and how much timsTOF its pretraining data contains.
- New papers since 10-30: the arXiv API rate-limited the sweep (HTTP 429), so one may be missing. The HF Hub showed
  no new MassSpecGym checkpoints.
- What the 0.540 leaders do. None has published.
