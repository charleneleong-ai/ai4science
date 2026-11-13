## Highlights
**Status 2026-11-13 (AEDT):** best **0.487** (extlib-fork), from 0.115 at the first submission. It is Hongbin Gu's
public 0.487 notebook *casmi26-v589-alex0484-extlib* run as our private kernel, unchanged except for a smoke commit
run: popwenman's Gold Apex V22 with an 80-epoch retrained spectrum model ensembled 0.5/0.5 with the public one, the
timsTOF mass prior, a documentation prior and 324,884 external library spectra (Takumu Hata's extlib plus harmonised
MassBank) for rank-1 promotion. Its +0.008 over the 0.479 dedup fork is one to two molecules of the ~124-molecule
public split, inside the ±0.016 between submissions, so it is the new base (PROMOTE), not a measured gain. 34
submissions scored; its kernel version 2, the same notebook with the learned isomer reranker on (`V8 = 'BOTH'`,
`RR8_ALPHA = 0.20`), scored 0.487 (no change, discarded) and version 3 (`RR8_ALPHA = 0.15`) is running. Two offline probes closed without a submission (MIST `msg_aug` encoder, an external
GNPS+MassSpecGym library). The leaderboard top is 0.540 (2026-11-12).

| | public MRR@25 | submission |
|---|---|---|
| current best (inside the noise of 0.479) | **0.487** | extlib-fork |
| dedup fork (inside the noise of 0.477) | 0.479 | dedup-fork |
| merged fork | 0.477 | merged-fork |
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
- **The public lineage still moves the line:** the 0.487 extlib notebook (two neural scorers, mass prior,
  documentation prior, 324,884 external library spectra) reproduced at 0.487 (extlib-fork), +0.008 over 0.479; each
  step is inside the noise, but the public line rose 0.477 → 0.487 in twelve days.
- **More copies of the same neural scorer hurt:** averaging the public spectrum model with three retrains of it
  (same recipe, new seeds) on ranks 2–25 lost 0.016 (seed-blend, 0.463). The seeds rank like the public model
  (Spearman 0.91–0.93) and two of three validate lower (0.481, 0.487, 0.465 against 0.484).
- **The public notebooks' mass prior moved nothing:** penalising candidates more than 3 ppm from the timsTOF
  mass-error centre on ranks 2–25, plus its rank-1 rule, scored 0.479 on the 0.479 base (mass-prior). The band fits
  the truth (99.3% of 2,138 timsTOF train answers inside it), so it is safe but does not reach enough answers.
- **The 0.479 rebase is not a measured gain:** +0.002 over 0.477 is below G1 (+0.005) and inside the source's own
  run-to-run spread (~0.002); it is kept because it is the current public line, not because it was shown to help.
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
| 2026-10-31 | merged-fork | the public 0.477 notebook (leader-safe neural re-rank + fragment check on rank 1), unchanged | 0.477 |
| 2026-11-01 | dedup-fork | the public 0.479 notebook (documentation prior + InChIKey14 dedup on ranks 2–25), unchanged; inside the noise | 0.479 |
| 2026-11-06 | seed-blend | dedup-fork + ranks 2–25 neural score averaged over the public scorer and three seed retrains | 0.463 |
| 2026-11-07 | mass-prior | dedup-fork + timsTOF precursor-mass prior on ranks 2–25 and a mass rank-1 rule (from Gold Apex V22) | 0.479 |
| 2026-11-12 | extlib-fork | the public 0.487 notebook (Gold Apex V22 + 80-epoch scorer ensemble + external spectral libraries), unchanged; inside the noise | **0.487** |
| 2026-11-13 | isomer-rerank | extlib-fork + learned isomer reranker (`RR8_ALPHA = 0.20`); no change, discarded | 0.487 |

## Task
[CASMI 2026](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra): rank
up to 25 structures per molecule from its LC-MS/MS spectra (400 molecules, 1213 spectra). Metric:
MRR@25 on the InChIKey connectivity block. Notebook-only, internet off; only the public split is
visible.

## Approach
The current best is an unchanged fork of the public notebook
[CASMI26 LB 0.479 | MIDX Exact Dedup Fill](https://www.kaggle.com/code/alexchilton/casmi26-lb-0-479-midx-exact-dedup-fill)
(alex chilton, running Baidalin Adilzhan's *midx exact dedup fill* unchanged),
[`notebooks/fork479/fork479.ipynb`](https://github.com/charleneleong-ai/ai4science/blob/feat/casmi-fork/casmi/notebooks/fork479/fork479.ipynb).
It is the 0.477 merged notebook
([`fork477`](https://github.com/charleneleong-ai/ai4science/blob/feat/casmi-fork/casmi/notebooks/fork477)) with a
final layer on ranks 2–25: a documentation prior (PubChem popularity 0.10, a spectral-library bonus 0.05 where the best
library match reaches 0.87, a natural-product bonus 0.02), the old leader pinned at rank 2 after a fragment-check move,
and deduplication by the metric's InChIKey14 key, with freed slots filled by scored outside candidates. Underneath, it
runs the pinned 0.465 pipeline once ([`fork464`](https://github.com/charleneleong-ai/ai4science/blob/feat/casmi-fork/casmi/notebooks/fork464)):
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
| documentation prior + leader pin + InChIKey14 dedup on ranks 2–25 (merged-fork → dedup-fork) | +0.002, inside the noise |
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
- Analogue propagation across test molecules on ranks 2–25 (analogue-fork, 2026-10-31): 0.477, no change from the
  0.477 base. Anchors with library match ≥ 0.87 boost candidates similar to their leader in molecules at spectral
  cosine ≥ 0.70 (at most 113 / 400 reachable). On labelled natural products such neighbours share structure 46% of
  the time; not enough to move the leaderboard at λ 0.10.
- Seed ensembles of the public spectrum model (seed-blend, 2026-11-06, 0.463 against 0.479): three same-recipe
  retrains on the A100 (16 × 500 steps, about 1 h 26 min each), blended with the public scorer by per-molecule
  z-score at equal weight. The run was healthy (all 400 molecules, 300 with every member scoring). Longer schedules
  had already scored lower in public (40 × 2000 steps: 0.453 against 0.472). Not retried with other weights.
- The timsTOF precursor-mass prior (mass-prior, 2026-11-07, 0.479 against 0.479): romanrozen's layer from popwenman's
  Gold Apex V22 with its constants, alone on our base. Ranks 2–25 lose up to 0.15 when a candidate's exact mass sits
  more than 3 ppm from −1.5 ppm of the neutral precursor mass; an in-band rank 2–3 replaces a leader more than 5 ppm
  off when the best library match is below 0.9. On 2,138 timsTOF train spectra the true structure sits at median
  −1.9 ppm, 99.3% inside the band; on the 12 smoke lists 59 of 288 rank 2–25 candidates fall outside it and the rank-1
  rule moved none. V22's other changes (scorer blend, prior weights, LOTUS backfill) were not tried.
- The extlib fork (extlib-fork, 2026-11-12, 0.487 against 0.479): Hongbin Gu's public notebook, which huichieh ran
  unchanged for 0.487 (submission 58185447), forked with one change, `COMMIT_MODE = 'smoke'`, so the gate sees every
  stage on the 12 visible molecules. The gate `gate_fork487.py` checked the OURS bundle and the ensemble line, the mass
  prior (394 / 400 lists touched), the documentation prior and the extlib stage (324,884 spectra loaded, 0 matches on
  the train-copy smoke molecules). Submission 58272211 reproduced 0.487.
- The isomer reranker (isomer-rerank, 2026-11-13, 0.487 against 0.487): kernel version 2 turns on the notebook's
  learned reranker (`V8 = 'BOTH'`, `RR8_ALPHA = 0.20`), which reorders same-formula isomers on lists whose library
  match is below 0.9. The gate `gate_fork487_v2.py` saw it run on all 400 lists, but the 12 smoke lists are library
  matches, so it reranked none and its reach on the hidden test was not visible. Submission 58292050 scored 0.487, no
  change: discarded. ozertuu's public "[0.500+]" notebook uses these exact settings; its 0.500 is a target, not a score.
- MIST `msg_aug` encoder for the de novo route (offline, 2026-11-01, no submission): mean fingerprint Tanimoto 0.385 on
  the 50 FRIGID panel molecules against a pre-registered gate of 0.40 (the old `mist_msg.pt`: 0.261). The `noaug`
  control is unusable (0.19 even on its own training spectra), so the source of the gain is not measured. The
  MIST/FRIGID line stays closed.
- An external GNPS + MassSpecGym library ([Zenodo 20042904](https://zenodo.org/records/20042904), offline, 2026-11-01,
  no submission): 3,096 of its 41,210 structures (13,814 spectra) are not already in train, and they match 0 / 400 test
  molecules at cosine ≥ 0.5 within 20 ppm in the same ion mode. Caveat: every visible test spectrum matches train at
  ≥ 0.87, so the visible test looks like train stand-ins, not the hidden test; for novel Enveda-180-like structures
  exact reach is about 0 by construction.
- PubChem promotion at rel > 100 on E006: not run; its public sibling already scored 0.444 against 0.446.
- FRIGID de novo (offline, 2026-10-18, no submission): 0 / 50 exact in the top 10 on Enveda-180 molecules with no
  public spectra, even with the true formula (bar fixed beforehand: 5 / 50). Its MIST encoder's fingerprints reach
  Tanimoto 0.26 to the truth on Enveda spectra (0.63 on MassBank); the decoder is not the bottleneck.
- FRIGID driven by the fold-safe FPNet instead: FPNet's fingerprints are no closer (0.245, own bit space) and use a
  different bit space; retraining a 4096-bit head (est. 10–15 T4 hours) was not justified.
- Excluding the train structure behind a copied spectrum, and similarity filters around it: the copies
  are a property of the visible test (train examples), not the hidden test the leaderboard scores.

## Next
The base is now the 0.487 extlib fork. Its kernel version 2 turned on the notebook's learned isomer reranker at a small
weight (`V8 = 'BOTH'`, `RR8_ALPHA = 0.20`; it reorders only same-formula isomers on lists whose library match is below
0.9) and scored 0.487, no change; version 3 lowers the weight to 0.15 and is running. The seed blend (0.463) closed the same-model ensemble route, and the Gold Apex V22 mass prior
alone scored 0.479 (no change). The daily watch of public notebooks runs until the deadline and flags a claim of 0.485
or more. One public notebook puts the public split at about 125 molecules, ±0.016 between submissions by chance
(unverified). Final picks (two, by 2026-12-14 23:59 UTC): the extlib fork (pick A) and a diverse line as pick B, the
0.479 dedup fork or the 0.477 merged fork, so the two picks do not share the external-library promotion. The ICEBERG atlas was reviewed on
2026-11-03 and skipped. The one resource left unused by
every public notebook is the [ICEBERG PubChem atlas](https://iceberg-ms.mit.edu/) (predicted spectra for 91.4M PubChem
structures), which the host allowed as a private dataset on 2026-10-28
([discussion 745134](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/745134));
free open pre-built datasets are allowed generally
([753662](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753662)). Analogue propagation scored 0.477 (no change). No offline panel can rank changes at this
level: sanpier's labelled panel ranks the 0.477 recipe below 0.465, so only the leaderboard measures, at about ±0.010
per comparison. Before the 0.477 fork, the base was the 0.465 data fork. Single settings of the E006 family look tuned (promotion threshold swept upstream,
engine-2 weight 0.8 lost 0.009, forward models already at 1.0, a neural re-rank at weight 0.50 lost 0.021), and the
de novo route failed offline. The public
lineage has gained +0.044 in a week (0.421 → 0.465), so a daily check of new public notebooks, forking any credible
improvement that reproduces, is the cheapest lever left. The leaderboard top is 0.540 (then 0.539 and 0.534, 2026-11-12). A fork costs about
2.5–3 GPU hours of the weekly 30 h quota, including the scoring rerun.
