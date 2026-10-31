# CASMI 2026 SOTA delta sweep: 2026-10-17 → 2026-10-30

This is a delta on the [2026-10-17 sweep](2026-10-17-casmi-sota-sweep.md). It covers only what is new or newly
relevant.

Raw notes:
- [`kaggle_discussion.md`](../research_notes/2026-10-30-casmi-sota-delta/kaggle_discussion.md): 36 new forum
  topics, the leaderboard and the public notebooks;
- [`methods.md`](../research_notes/2026-10-30-casmi-sota-delta/methods.md): arXiv, Hugging Face, GitHub and
  Zenodo.

Paper numbers are the authors' own. Forum and notebook numbers are participants' single public-LB submissions and
are not verified here. Our position: a fork of the public 0.477 notebook ([`c467db3`](https://github.com/charleneleong-ai/ai4science/commit/c467db3)),
expected to score about 0.477 on the public LB.

## Leaderboard and field

- **The top moved from 0.495 to 0.517 in 13 days.** On 2026-10-30 the top five were:

  | Team | Score |
  |---|---|
  | Oleksandr Parfionov | 0.517 |
  | Savko | 0.513 |
  | Shehab Alshehabi | 0.501 |
  | pikachu | 0.499 |
  | CHU1PC | 0.495 |

  The 0.040 gap from 0.477 is about 5 public molecules.
- **The public field is one pipeline.** A census in
  [kunaldesale2408's notebook](https://www.kaggle.com/code/kunaldesale2408/enveda-casmi-2026-mid) (2026-10-27):
  - 103 of the 112 public notebooks at 0.460 or above attach dmitriigluzdov's v4g assets;
  - 130 teams sit at exactly 0.472, 169 at 0.468 and 190 at 0.464;
  - the gold line was 0.477 (rank 18);
  - "No member of the top ten teams is the author of any of the 974 public notebooks."

  None of the leaders has posted a method.
- **Public-LB noise is now quantified**
  ([753511](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753511)):
  - With about 125 public molecules and 10% of molecules changing, a change with no real effect still swings by
    0.010 or more in 47% of comparisons.
  - Picking the best of 7 such variants gains about +0.019 on the public LB and nothing on the private LB.
  - A paired offline comparison on 600 molecules has a standard error of about 0.0065.

## Host rulings since 2026-10-17

| Ruling | Thread | Consequence |
|---|---|---|
| Automated, unsupervised processing **across test spectra** at inference is allowed, e.g. "using a confidently identified test molecule to help rank candidates for similar test molecules" | [752297](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/752297) (Marie Killian, 10-26), [751373](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/751373) (David Healey) | Opens transductive analog propagation, which no public notebook has documented (**unverified**) |
| CC BY-NC open data allowed generally; NIST and its derivatives still banned | [750330](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/750330) | HMDB-trained models and similar are fine |
| MS2DeepScore 2.0 and FIORA-OS allowed; SIRIUS/CSI:FingerID **not** allowed; PubChem-only chemical language models allowed as priors; the coleygroup MIST checkpoints (Zenodo 23303289) allowed | 751373 | New encoder option (see Methods) |
| Intermediates need their generating code plus the data version or date, not bit-for-bit reproduction | 751373 | |
| "We reserve the right to remove submissions that are not reproducible from the leaderboard, regardless of ranking" | [752541](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/752541) | **Risk for our fork** (below) |
| **Still unanswered:** whether participant-published checkpoints and intermediates without training code are compliant. Named assets: megayak simulated ranker rows, prvsiyan fingerprint caches, ahmedberatozer PubChem tier, dmitriigluzdov NP knowledge table | [753775](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753775), [753662](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753662), [753900](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753900) | Our 0.477 fork depends on several of these, and on Saner Turhaner's neural scorer |

## Measured results from the public lineage

These are single public-LB submissions with an SE of about 0.014, so read them as direction only.

- **The 0.477 recipe itself.** busyaprime's submission table, reproduced in kunaldesale2408's notebook:
  - Fusion rule: 1/(3+b) + w/(3+n), with w = 0.25 at rank 1 and 0.75 at ranks 2-25. Here b is the pipeline rank
    and n the neural rank.
  - Adding the fragment check's leader moves from the 0.50-order list gave **0.477**, against 0.472 without them.
  - Swapping in Ahmed Berat Özer's longer-trained neural scorer gave **0.453 / 0.461**.
  - An "NX isomer reranker" at weight 2 gave **0.414**.
  - So a stronger isomer or neural model has *lost* every time it got more say.
- **dmitriigluzdov**
  ([Data Behind](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-the-data-behind-spectra-to-structures) §6):
  - removing PubChem-only candidates: 0.464 → 0.409;
  - 6 → 20 analogue parents for generation: 0.464 → 0.456;
  - five fold-safe FPNets instead of one: no change;
  - his own competition-trained neural re-order with leader-safe fusion: 0.459 against 0.465.
- **Class 3 is still not moving on the LB:**
  - [Kino Lutor](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/750241):
    on 130 molecules missing from the databases, the truth is rank 1 only 7 times, although the analogue generator
    *can* produce 100 of the 130. More parents still lost on the LB (above).
  - [blaise yonga](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753400):
    65 biotransformation rules generate the truth for 46% / 32% / 16.5% of molecules when the nearest analogue is
    capped at no limit / Tanimoto 0.70 / 0.55. Class-3 MRR is still only 0.080 / 0.056 / 0.030. Expected net gain
    +0.004 to +0.010, so not used.
  - [Bull and Tra Crew](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/751300):
    weeks of class-3 work gave flat or +0.003-0.004 on the LB.
- **Where the losses are.**
  [Pathik Patel's failure ledger](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/752735)
  (n=500, own honest split):

  | Bucket | Share |
  |---|---|
  | In no database | 21.6% |
  | In the mass window but outside the top 25 | 14.8% |
  | Ranked 2-25 | 39.4% |
  | Ranked 1 | 23.8% |

  Pool expansion to 97M structures gave +0.113 locally but +0.003 on the LB.

## Methods

| Item | What / primary link | Author-reported numbers | Offline / T4 | Reach on our pipeline | Effort |
|---|---|---|---|---|---|
| **FRIGID / MIST leakage correction** | The FRIGID MSG MIST encoder had been trained on synthetic spectra of PubChem molecules matched to test formulae. Retrained in [FRIGID README](https://github.com/coleygroup/FRIGID) and [arXiv 2604.16648 v3](https://arxiv.org/abs/2604.16648) (10-25) | Top-1 fell from 25.03% to **14.16%** (NPLIB1) and from 18.29% to **9.05%** (MSG) | n/a | Our FRIGID 0/50 is less surprising: the published bar was about half as high. De novo stays closed | none |
| **MIST checkpoints, including `msg_aug`** | [coleygroup/MIST](https://github.com/coleygroup/MIST) (MIT, training code) and [Zenodo 23303289](https://zenodo.org/records/23303289) (10-20). The aug model is MSG plus GLACIER-simulated spectra of 4M molecules (78M rows). Host-allowed | No timsTOF numbers | About 290 MB checkpoint, fine on T4 | Our closed MIST test used FRIGID's old `mist_msg.pt` (Tanimoto 0.26 on Enveda-180). Our own experiment card names the **encoder** as FRIGID's bottleneck: fed the true fingerprint, the decoder recovered 49/50 (inflated by pretraining) | 0.5 day offline probe |
| **Separability Is Not Alignability** | [arXiv 2610.37979](https://arxiv.org/abs/2610.37979) (10-27), [code](https://github.com/autonlab/msms-retrieval), no weights | MSG v1.5 MCES split, decontaminated: frozen DreaMS/MoLFormer alignment caps at MRR 0.227; pretraining on 5.9M GLACIER + FraGNNet simulated pairs plus MaxSim reaches **MRR 0.566, R@1 0.425**. Training on candidate negatives learns a spectrum-blind "library vs decoy" prior worth ΔMRR 0.12-0.14 | Training 4×A100; inference fits a T4 | Recipe for a much stronger re-ranker than Saner Turhaner's, but the busyaprime and dmitriigluzdov rows show stronger neural scorers *losing* on the LB | 3-5 days; too big now |
| **GEMS** | [arXiv 2610.23617](https://arxiv.org/abs/2610.23617) (10-16): learned bond-edit search over same-formula isomers, seeded from PubChem | MSG hit@1/10 30.2/41.3%; NPLIB1 29.6/50.8%. Uses the pre-fix MIST, so likely inflated | No code | Would edit near-miss pool members (the 14.8% + class-3 buckets) | blocked |
| **SPARC** | [arXiv 2609.35649 v2](https://arxiv.org/abs/2609.35649) (10-18): reference-guided test-time recalibration of a forward model for acquisition shift | No CASMI-like numbers | No code | Idea for timsTOF vs MSG-trained GLACIER/ICEBERG | 1-2 days, idea only |
| **Selective prediction for retrieval** (Aug, missed last time) | [arXiv 2603.10950](https://arxiv.org/abs/2603.10950) | Ranking-level margin and entropy best predict whether rank 1 is correct; fingerprint-bit uncertainty does poorly | trivial | A per-molecule gate for how much say the neural or fragment layer gets | 0.5 day |
| MassSpecGym leaderboard, ms-pred, DreaMS, HF Hub | [MassSpecGym](https://github.com/pluskal-lab/MassSpecGym) has no commits to main since 2026-05-08; ms-pred had one data-loader commit; HF has no new relevant checkpoints | none | | no change | |

## What changed since 2026-10-17

1. **Leaderboard:**
   - the top is now 0.517, and nobody in the top ten publishes;
   - the public lineage saturates at 0.477;
   - noise is now quantified: a single public comparison cannot resolve ±0.010.
2. **New legal lever:** transductive use of the test set (cross-molecule spectral similarity) is explicitly
   allowed as of 2026-10-26.
3. **New legal model:** the MIST checkpoints (MIT, with training code), including the GLACIER-augmented `msg_aug`.
4. **De novo's published numbers halved** after the leakage fix (FRIGID MSG 18.3% → 9.1%). With the forum's
   class-3 negatives, this confirms class-3 generation is not a near-term lever.
5. **"More model" keeps losing on the LB:**
   - stronger neural scorers lose 0.011-0.019;
   - an isomer reranker loses 0.058;
   - more generation parents lose 0.008.

   The gains that held were small, conservative layers on ranks 2-25 plus rare, gated rank-1 moves.
6. **Eligibility risk is now concrete.** Non-reproducible submissions can be removed from the LB, and the status
   of the shared public assets is unanswered.

## Does it change the plan

Yes, in two ways:
- **Stop scaling models and start using the test set's own structure.** Transductive propagation is the only
  newly opened lever that no public notebook documents. It acts on the 39% "ranked 2-25" bucket and on analogue
  families, which "real or potential natural molecules and their analogs" suggests are common in this test.
- **Plan final selection for eligibility.** Keep one of the two final submissions on a pipeline whose every weight
  and intermediate has public generating code, until the host answers 753775.

Ranked candidate experiments, one change each. Each is gated offline before it costs a submission, and a
submission counts as a gain only if it clears about +0.010.

| # | Experiment (one change) | Why now | Zero-submission gate | Expected public-LB gain | Effort |
|---|---|---|---|---|---|
| 1 | **Transductive analogue propagation across test molecules.** Anchors are test molecules whose rank 1 is confident (library cosine ≥ 0.87 or a large fusion margin). For every other test molecule, add a ranks 2-25 bonus of max over anchors of `spec_sim(test_i, anchor) × Tanimoto(candidate, anchor_rank1)`. Rank 1 is protected, with the same leader-safe form as the 0.477 rule | Host-approved 10-26 (752297); no public notebook documents it; the test is natural products and analogues; acts on the largest loss bucket (ranked 2-25) | On Enveda-180 molecules in train without public spectra (Aryan's class-2-like split, identity-disjoint folds), measure: (a) the share of molecules with a within-set neighbour at spectral cosine ≥ 0.7; (b) the paired MRR delta with SE. Submit only if (a) ≥ 15% and the delta exceeds 2 × SE | +0.004-0.015 (speculative; depends on how clustered the hidden test is) | 1-1.5 days |
| 2 | **Confidence-gated neural weight.** Use w_mid = 0.75 only where the pipeline's top-1 vs top-2 fused-score margin is below a threshold fixed offline, else 0.50. Nothing else changes | The 0.50 vs 0.75 choice moved the LB by about 0.005 and our own 0.50 run scored 0.444. Ranking-margin confidence is the best per-query correctness signal ([2603.10950](https://arxiv.org/abs/2603.10950)) | On the same panel, fit the threshold on one fold, then report the paired delta on the others. Submit only if it exceeds 2 × SE (≈ 0.013) | +0.000-0.008 | 0.5 day |
| 3 | **Re-probe the MIST `msg_aug` encoder** ([Zenodo 23303289](https://zenodo.org/records/23303289)) on the 50 Enveda-180 FRIGID panel molecules, as a direct replacement for the old `mist_msg.pt` | Our closed test blamed the encoder (Tanimoto 0.26), and the new checkpoint is trained on 78M GLACIER-simulated spectra with the leak fixed | Offline only: predicted-fingerprint Tanimoto on the 50. If it is at least 0.40, reopen the gated FRIGID insert at ranks 3-5 (2026-10-17 plan #1) and as an FPNet A/B. Otherwise MIST stays closed for good | 0 unless the gate passes; then +0.005-0.010 (class 3, speculative) | 0.5 day probe |

Not ranked:
- **Separability-style retraining** (5.9M simulated pairs): too big and against the LB evidence.
- **GEMS:** no code.
- **SPARC:** no code.
- **More analogue parents or biotransformation generation:** measured negative.

**Run alongside these, at no submission cost:**
- an audit of which attached assets in the 0.477 fork lack public generating code, so a reproducible fallback
  exists for final selection;
- a spectrum-blind control of the engine-2 ranker, which tells us how much of its lift is the decoy prior that
  2610.37979 describes.

## Not established

- The class shares, which the host will not release, and how clustered the hidden test is. The latter sets
  experiment 1's reach, and only the hidden run can measure it.
- Any leader's method.
- Whether shared public assets without code will be ruled ineligible.
- MIST `msg_aug` performance on timsTOF spectra.
- All T4 runtimes for new items.
