## Highlights
**Status 2026-10-31:** best **0.477** (merged-fork), from 0.115 at the first submission. It is Busya PRIME's public
0.477 notebook run unchanged as our private kernel: the 0.465 data-fork pipeline, then a leader-safe neural re-rank
(spectrum model at weight 0.25 for rank 1, 0.75 for ranks 2–25) and Roman Tamrazov's fragment check on rank 1.
28 submissions scored; fork477 v3 (analogue propagation across test molecules) is scoring. The leaderboard top is 0.517.

| | public MRR@25 | submission |
|---|---|---|
| current best | **0.477** | merged-fork |
| data fork | 0.465 | data-fork |
| PubChemLite fork | 0.456 | pubchemlite-fork |
| E006 fork | 0.446 | e006-fork |
| fork of the 0.421 notebook | 0.415 | fork-v1 |
| our own pipeline's best | 0.356 | v19 |
| first submission | 0.115 | v1 |

- **Rebasing on the public notebooks was the biggest gain:** +0.059 for the 0.421 lineage, then +0.031 more for E006
  (engine-2 rank fusion, ChEBI/LipidMaps candidates, PubChem promotion at rel > 150). Its library search over train,
  analog propagation and LightGBM ranker use train matches as evidence, where our pipeline excluded them. Two later public forks added +0.019 (PubChemLite in engine 2, then proxy re-ordering and fold-safe FPNets).
- **A neural re-rank helps only when it cannot move rank 1:** at weight 0.25 for rank 1 and 0.75 below it, plus a
  fragment check on rank 1, it added +0.012 (merged-fork); the same model at a flat 0.50 lost 0.021 (neural-fork).
- **A learned spectrum model was the biggest lever:** FPNet in place of fragment matching, +0.115.
- **Forward models were the second:** ICEBERG then GLACIER, +0.087 from v11 to v15, most of it from
  GLACIER (+0.028) and from weighting both at 1.0 instead of 0.5 (+0.029).
- **Candidate-pool and prior tuning moved the score by under 0.01**, as did removing the copied
  train structures, even though that changed 267 top-1 picks.
- **No local estimate exists.** The first local estimate (0.664) scored 0.115; every number here is
  a leaderboard score. Train has no test-like subset to calibrate one on: 41 isomer pairs in 2.54M
  rows share a spectrum, and 7 split across the fold-safe FPNet's training the way the test's pairs do.

| date | version | what changed | public MRR@25 |
|---|---|---|---|
| 2026-10-07 | v1 | fragment explainability over train | 0.115 |
| 2026-10-08 | v6 | FPNet | 0.250 |
| 2026-10-09 | v12 | + ICEBERG | 0.286 |
| 2026-10-10 | v14 | + GLACIER | 0.320 |
| 2026-10-10 | v15 | forward-model weights 1.0 | 0.349 |
| 2026-10-12 | v19 | + second PubChem pass, skipping candidates similar to the copied structure | 0.356 |
| 2026-10-17 | fork-v1 | the public 0.421 notebook, unchanged, as our kernel | 0.415 |
| 2026-10-17 | e006-fork | the public 0.446 E006 notebook, unchanged, as our kernel | 0.446 |
| 2026-10-20 | pubchemlite-fork | the public 0.456 notebook (+ PubChemLite in engine 2), unchanged | 0.456 |
| 2026-10-24 | data-fork | the public 0.464 notebook (v26: proxy re-ordering, fold-safe FPNets, candidate-aware promotion), unchanged | 0.465 |
| 2026-10-31 | merged-fork | the public 0.477 notebook (leader-safe neural re-rank + fragment check on rank 1), unchanged | **0.477** |

## Task
[CASMI 2026](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra): rank
up to 25 structures per molecule from its LC-MS/MS spectra (400 molecules, 1213 spectra). Metric:
MRR@25 on the InChIKey connectivity block. Notebook-only, internet off; only the public split is
visible.

## Approach
The current best is an unchanged fork of the public notebook
[LB 0.477 · Two 0.472 Fixes, Merged the Right Way](https://www.kaggle.com/code/busyaprime/lb-0-477-two-0-472-fixes-merged-the-right-way),
[`notebooks/fork477/merged_fork.ipynb`](https://github.com/charleneleong-ai/ai4science/blob/feat/casmi-fork/casmi/notebooks/fork477/merged_fork.ipynb).
It runs the pinned 0.465 pipeline once ([`fork464`](https://github.com/charleneleong-ai/ai4science/blob/feat/casmi-fork/casmi/notebooks/fork464)):
library search of every test spectrum against train with analog propagation, a public LightGBM ranker, PubChem,
PubChemLite, ChEBI and LipidMaps candidates, a second engine's ranks fused in, ICEBERG and GLACIER re-scoring with
proxies, fold-safe FPNets and a candidate-aware rank-1 promotion. On top it re-ranks the final 25 with a public
spectrum model by reciprocal-rank fusion, at weight 0.25 for rank 1 and 0.75 for ranks 2–25, and lets an untrained
fragment check (Spectral Parallax) replace rank 1, its moves taken from the 0.50 order. Our own pipeline,
[`kaggle_submission.py`](https://github.com/charleneleong-ai/ai4science/blob/main/casmi/notebooks/kaggle_submission.py)
(best 0.356), is the record of what we built before the forks. Inputs and checkpoints:
[experiment card](https://github.com/charleneleong-ai/ai4science/blob/main/casmi/docs/experiment-card.md);
summary: [report](https://github.com/charleneleong-ai/ai4science/blob/main/casmi/docs/report.pdf).

## Results
Generated by `casmi wandb-sync`.

## Full run log
Generated by `casmi wandb-sync`.

## What we learned
| lever | measured gain |
|---|---|
| the public notebook's library search + ranker in place of our pipeline (v19 → fork-v1) | +0.059 |
| E006 lineage in place of the 0.421 notebook (fork-v1 → e006-fork) | +0.031 |
| PubChemLite in engine 2 + N1A + popularity 0.15 (e006-fork → pubchemlite-fork) | +0.010 |
| proxy re-ordering, fold-safe FPNets, candidate-aware promotion (pubchemlite-fork → data-fork) | +0.009 |
| leader-safe neural re-rank + fragment check on rank 1 (data-fork → merged-fork) | +0.012 |
| FPNet in place of fragment explainability | +0.115 |
| popularity shortlist instead of nearest-mass cut | +0.086 |
| GLACIER beside ICEBERG | +0.028 |
| forward-model weights 1.0 instead of 0.5 | +0.029 |
| ICEBERG (on CPU, partial coverage) | +0.024 |
| second PubChem pass with a copied-structure similarity filter (v15 → v19; threshold post hoc) | +0.007 |

**Closed — do not retry:**
- Nearest-mass candidate cut: arbitrary among same-formula isomers (0.049).
- Forward-model weights 1.5: no gain over 1.0 (0.346). ICEBERG at 0.5 beside GLACIER 1.0: 0.333 (v21).
- Spectrum-driven second PubChem pass without a filter (v18, 0.353): +0.164 MRR@25 on held-out
  train molecules, +0.004 on the leaderboard. The production FPNet memorised the copied spectra and
  pulled the pass toward the copied structures (median Tanimoto of changed top-1 picks 0.26 → 0.42).
  Skipping candidates with Tanimoto ≥ 0.5 to the copied structure (v19, 894 dropped) took it to 0.356.
- Choosing the second pass with the fold-safe FPNet: it trained on the copied rows of 317 / 400 test
  molecules, so it memorised them too; not submitted.
- GLACIER on `[M+H]+` only: within noise (0.345).
- Surrogate spectra for the forward models: reach 2 / 400 more molecules, below noise; not submitted.
- LightGBM ranker over our own features (FPNet z / rank, popularity, natural-product flags, patents,
  ppm error, isomer count): ties the hand-set score in 5-fold grouped CV on 1,881 held-out train
  molecules (0.296 vs 0.296, Δ +0.001 [−0.007, +0.008]); not submitted.
- PubChem channel over the top 10,000 pass-1 rows instead of 5000 (fork-v2): 0.415, no change.
- Engine-2 weight in the e006 fork's rank fusion 0.8 instead of 0.6 (e006-v3): 0.437, −0.009.
- The public 0.451 notebook at its current version 53 (spectra-fork): 0.445. Kaggle serves only the current
  version; upstream's scored version 52 (0.451) is not available to fork.
- A neural spectrum-model re-rank of the final 25 at weight 0.50 (neural-fork): 0.444, −0.022 against the same
  notebook's 0.466 at weight 0.25. The re-rank adds at most +0.001 to the 0.465 base; more weight hurts.
- PubChem promotion at rel > 100 on E006: not run; its public sibling already scored 0.444 against 0.446.
- FRIGID de novo (offline, 2026-10-18, no submission): 0 / 50 exact in the top 10 on Enveda-180 molecules with no
  public spectra, even with the true formula (bar fixed beforehand: 5 / 50). Its MIST encoder's fingerprints reach
  Tanimoto 0.26 to the truth on Enveda spectra (0.63 on MassBank); the decoder is not the bottleneck.
- FRIGID driven by the fold-safe FPNet instead: FPNet's fingerprints are no closer (0.245, own bit space) and use a
  different bit space; retraining a 4096-bit head (est. 10–15 T4 hours) was not justified.
- Excluding the train structure behind a copied spectrum, and similarity filters around it: the copies
  are a property of the visible test (train examples), not the hidden test the leaderboard scores.

## Next
The base is now the 0.477 merged fork. Running: fork477 v3, which adds analogue propagation across test molecules
(allowed by the host on 2026-10-26; at most 113 / 400 molecules reachable; no offline panel can score it, since
sanpier's labelled panel ranks the 0.477 recipe below 0.465). It counts only if it beats 0.477 by more than ~0.010.
Before that, the base was the 0.465 data fork. Single settings of the E006 family look tuned (promotion threshold swept upstream,
engine-2 weight 0.8 lost 0.009, forward models already at 1.0, a neural re-rank at weight 0.50 lost 0.021), and the
de novo route failed offline. The public
lineage has gained +0.044 in a week (0.421 → 0.465), so a daily check of new public notebooks, forking any credible
improvement that reproduces, is the cheapest lever left. The leaderboard top is 0.504. The weekly GPU quota (30 h)
resets 2026-10-31 00:00 UTC; a fork costs about 2.5–3 GPU hours including the scoring rerun.
