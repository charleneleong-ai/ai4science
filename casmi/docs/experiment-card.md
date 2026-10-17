# Experiment card — CASMI 2026 molecule identification

Pattern: [*A generated protein is not a result*](https://rewirebio.io/blog/a-generated-protein-is-not-a-result/),
applied to a retrieve-and-rank campaign. Status as of 2026-10-17: the current best (0.446) is an
unchanged fork of the public E006 notebook, [`notebooks/fork446/e006_fork.ipynb`](../notebooks/fork446/e006_fork.ipynb);
the previous best (0.415) is an unchanged fork of the 0.421 notebook,
[`notebooks/fork421/v4n_fork.ipynb`](../notebooks/fork421/v4n_fork.ipynb); our own pipeline (best 0.356) is [`notebooks/kaggle_submission.py`](../notebooks/kaggle_submission.py). History:
[PR #80](https://github.com/charleneleong-ai/ai4science/pull/80); every scored row:
[`verdicts.yaml`](experiments/casmi/verdicts.yaml).

## Target — frozen

| | |
|---|---|
| task | rank up to 25 candidate structures (SMILES) per molecule from its LC-MS/MS spectra |
| metric | MRR@25 on the tautomer-canonical InChIKey connectivity block (first 14 characters) |
| test set | `test.parquet`: 400 molecules, 1213 spectra; answers hidden |
| scored on | the **public** leaderboard split only; the private split is not visible until the competition closes |
| rules | code competition: notebook only, internet off; about 5 submissions a day |

## Current best — e006-fork, public MRR@25 0.446

| | |
|---|---|
| kernel | `charyeezy/casmi-e006-fork` (private), kernel version 2, submission 57298936; reproduces the source's public 0.446 exactly |
| source | [CASMI26 \| E006 + PubChem Promo rel>150 \| LB 0.446](https://www.kaggle.com/code/seyitkaangunes/casmi26-e006-pubchem-promo-rel-150-lb-0-446), byte-identical as of 2026-10-17 (md5 `9829a9667030e7044027f0cd33fd7f5d`); only [`kernel-metadata.json`](../notebooks/fork446/kernel-metadata.json) differs (our id, private, the GPU image) |
| lineage | flexonafft v80 recipe (public 0.436) → E006: engine-2 rank fusion with the `rank2` hedge, ChEBI/LipidMaps candidates, simulated ranker rows (0.442) → PubChem promotion's structure-evidence clause at `rel > 150` instead of 250 (0.446) |
| ranker | engine 1: `ranker_v4b.pkl`, public LightGBM ranker, 160 features, 4 boosters; engine 2: two gradient-boosted rankers fitted in the notebook; none retrained by us |
| datasets | 16 attached unchanged, as listed in [`kernel-metadata.json`](../notebooks/fork446/kernel-metadata.json): `prvsiyan/casmi26-ranker-features`, `prvsiyan/casmi26-fp-models-v2`, `dmitriigluzdov/casmi26-pubchem-popularity-prior`, `megayak/casmi26-simulated-ranker-rows`, `ahmedberatozer/casmi26-{fpnet-full1, glacier, iceberg, pubchem-tier, v2-pool, v3-models, v4b-models}`, `prvsiyan/chebi-lipidmaps-casmi26`, `prvsiyan/coconut-casmi26-candidates`, `hengck23/hengck23-dreams-enveda-casmi26`, `metric/rdkit-2026-3-3-wheel`, `eveliaveldrine/rdkit-wheels` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461` (the one fork-v1 uses), NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_fork446.py` (outside the repo) passed on the commit run (1183 s): engine ready; ranker 160 features / 4 boosters; base errors 0; ICEBERG scored 11 / 12 and GLACIER 11 / 11 smoke molecules; engine-2 cell applied on 400 rows; output 400 rows. Hidden-test scoring took about 4 h |
| failed attempt | kernel version 1 kept the source's metadata, which names Kaggle's CPU image (`gcr.io/kaggle-images/python`): the notebook's preflight failed (`nvidia-smi` exit 12, no GPU) and it was not submitted |

The settings were fixed by the notebook's authors before we ran it; we chose none of them. Its
author selected `rel > 150` on the public leaderboard (a sweep on two base pipelines), so 0.446 may be
optimistic on the private split. Next: one setting per submission against 0.446 (a promotion
threshold below 150; the forward-model weights), with about 26 of 30 weekly GPU hours left until
2026-10-24 00:00 UTC.

## Prior best — fork-v1, public MRR@25 0.415

| | |
|---|---|
| kernel | `charyeezy/casmi-v4n-fork` (private), kernel version 1, submission 57288940 |
| source | [CASMI26 \| v4n Fusion + PubChem on public 0.421](https://www.kaggle.com/code/huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421) (public score 0.421), byte-identical as of 2026-10-13; only [`kernel-metadata.json`](../notebooks/fork421/kernel-metadata.json) differs (our id, private) |
| ranker | `ranker_v4b.pkl`: public LightGBM ranker, 160 features, 4 boosters; not retrained |
| datasets | attached unchanged: `prvsiyan/casmi26-fp-models-v2`, `dmitriigluzdov/casmi26-fold-safe-fpnet`, `dmitriigluzdov/casmi26-pubchem-popularity-prior`, `prvsiyan/casmi26-ranker-features`, `megayak/casmi26-simulated-ranker-rows`, `ahmedberatozer/casmi26-{fpnet-full1, glacier, iceberg, pubchem-tier, v2-pool, v3-models, v4b-models}`, `prvsiyan/chebi-lipidmaps-casmi26`, `prvsiyan/coconut-casmi26-candidates`, `metric/rdkit-2026-3-3-wheel` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461`, NVIDIA T4, internet off |
| run evidence | commit run 1286 s; ranker loaded; base lists 0 errors; PubChem join injected 600 rows in 12 smoke molecules, 0 failures; GLACIER scored 11 / 12 smoke molecules |

The fork's settings were fixed by its author before we saw any result; we chose none of them after
scoring. Its one change so far, `PC_N1=10000` (fork-v2), was chosen from our recall probe (255 against
236 / 280 held-out answers) before it was submitted, and scored 0.415: no change, discarded.

## Inputs to our own pipeline

| input | identity | used for |
|---|---|---|
| `train.parquet` (competition) | 275,798 unique structures | candidates; copied-spectrum detection |
| [`casmi26-pubchem-tier`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-pubchem-tier) | 105,875,489 structures; updated 2026-09-23 | candidates |
| [`casmi26-pubchem-popularity-prior`](https://www.kaggle.com/datasets/dmitriigluzdov/casmi26-pubchem-popularity-prior) | row-aligned with the tier; PubChem Extras 2026-09-27; updated 2026-10-05 | shortlist and prior |
| [`casmi26-v2-pool`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-v2-pool) | `fp_bits.npy`, 10,226 selected bits; updated 2026-09-21 | FPNet fingerprints |

**Not pinned:** Kaggle attaches each dataset at its latest version. The dates above identify what
the scored runs used; a dataset updated later would change a re-run silently.

## Models and tools — exact

| component | checkpoint / version | source |
|---|---|---|
| FPNet (spectrum → fingerprint) | `fpnet_full1.pt` + `fpnet_0.pt` (FPNet A), logits averaged | [`casmi26-fpnet-full1`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-fpnet-full1) (2026-09-29), [`casmi26-v4b-models`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-v4b-models) (2026-09-26) |
| ICEBERG (forward model) | ICEBERG 2.1, MassSpecGym `msg_all`: `ckpt/gen.pt`, `ckpt/inten.pt` | [`casmi26-iceberg`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-iceberg) (2026-09-26) |
| GLACIER (forward model) | `ckpt/glacier.pt`, MassSpecGym checkpoint | [`casmi26-glacier`](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-glacier) (2026-09-28) |
| RDKit | 2026.3.3 in the notebook; 2025.03.6 inside the forward-model subprocesses (their requirement) | [`metric/rdkit-2026-3-3-wheel`](https://www.kaggle.com/datasets/metric/rdkit-2026-3-3-wheel); ICEBERG's wheels |
| runtime | Python 3.12, image `gcr.io/kaggle-private-byod/python@sha256:37c64f7d…a461`, NVIDIA T4 | [`kernel-metadata.json`](../notebooks/kernel-metadata.json) |

FPNet full1 was **trained on all of train**, including the train rows the test spectra copy (see
below). No model here was trained or fine-tuned by this work; every weight is a public checkpoint.

## Budget — every attempt in the denominator

- **22 submissions, 22 scored, 0 rejected.** 19 from our own pipeline, one an identical re-run of v19
  to measure noise (0.353 vs 0.356); every one is in the
  [score table](https://github.com/charleneleong-ai/ai4science/pull/80), including those not kept
  (0.049, 0.254, 0.345, 0.346), plus v18 (0.353, within noise), v19 (0.356, kept), v20 (0.354) and v21 (0.333) in
  [PR #85](https://github.com/charleneleong-ai/ai4science/pull/85). Two from the 0.421 fork: fork-v1
  (0.415, 57288940) and fork-v2 (0.415, 57290631, discarded). One from the E006 fork: e006-fork
  (0.446, 57298936, kernel version 2).
- **Two measurement kernels, not submissions:** the ranker training set (2000 held-out train
  molecules, 3 h 37 min) and the second-pass recall probe (300 molecules, 35 min; v1 replaced before
  it finished because an empty mass window would have crashed it).
- **18 notebook versions of our own pipeline run on Kaggle,** plus 2 of the 0.421 fork and 2 of the E006 fork; of our own pipeline's, the 3 not submitted were setup or reproducibility runs
  (a missing RDKit wheel, a push-from-repo check). The E006 fork's version 1 ran on the CPU image
  its source's metadata names, failed the notebook's GPU preflight and was not submitted.
- **One run was cut short by budget, not failure:** v12's ICEBERG on CPU scored 269 / 371
  coverable molecules in its 6 h cap after the weekly GPU quota ran out.

## Gates and filters — and when each was fixed

| gate / filter | fixed | how decided |
|---|---|---|
| submit only if the run's log shows the models loaded, scores came back non-null, and the output is 400 × 25 | before v4 | process rule; no submission failed it |
| 5 ppm mass window; top 25 output | before v1 | metric and instrument convention |
| popularity shortlist of 1000 (from 300) | **after** scoring both (0.258 vs 0.261) | public leaderboard |
| drop the train structure behind each copied spectrum | before v11 was scored | analysis of the 0.421 notebook's picks, not of test labels |
| forward-model weights 1.0 (from 0.5; 1.5 tried) | **after** scoring all three | public leaderboard |
| GLACIER on all covered spectra (`[M+H]+`-only tried) | **after** scoring both | public leaderboard |
| second pass skips candidates with Tanimoto ≥ 0.5 to the copied structure | **after** v18 scored | v18's changed top-1 picks drifted toward the copied structures; 0.4 tried after (v20, 0.354) |
| the fork's settings (ranker, similarity 0.7 gate, PubChem channel, forward-model weights) | by its author, before we ran it | not ours; we cannot see what evidence they were tuned on |
| fork `PC_N1=10000` | **before** fork-v2 was submitted | our held-out recall probe |
| the E006 fork's settings (engine-2 fusion weights, `rank2` hedge, promotion at `rel > 150`) | by its authors, before we ran it | not ours; the threshold was swept on the public leaderboard |

Choices marked **after** were selected on the public leaderboard, so our pipeline's 0.356 is an
optimistic estimate of the private score; the forks' 0.415 and 0.446 involve none of ours, though the E006 threshold was tuned on the public split by its author. Differences under about ±0.005 — the run-to-run noise
of identical runs — were treated as no change.

## What is not measured

- **The private score.** Every number here is public-split.
- **Why fork-v1 scores 0.415 on our account against 0.421 public.** Both our runs agree, so the
  0.006 is not run-to-run noise; the E006 fork reproduced its public 0.446 exactly on the same account
  and day, so it is not a general account or day effect. The cause is not established.
- **What in E006 earns its +0.031 over fork-v1.** We ran it unchanged; the split between engine-2
  fusion, ChEBI/LipidMaps candidates, simulated ranker rows, the overlays and the threshold is unmeasured by us.
- **A per-molecule breakdown of the forks' gains.** The hidden test is not visible.
- **Any local estimate.** Train-derived labels were wrong (a local 0.664 scored 0.115), and
  public notebooks cannot supply labels: their downloadable outputs score only ~12 smoke-test
  molecules and fill the rest with a placeholder (`CCO` / `CN`), because full predictions run only
  inside Kaggle's scoring rerun (`IS_RERUN` in the 0.421 notebook). The earlier "three lineages
  agreeing on 12 / 400" reading was those placeholders. A test-like held-out set cannot be cut from train either: only 41 isomer pairs
  in 2.54M train rows share a spectrum (≥ 3 peaks), and 7 have one structure held out of the fold-safe
  FPNet and the other in its training, as test answers and copied structures are. The test's 1189 / 1213
  copied spectra were built deliberately, not by library duplication.
- **Why the data is shaped as it is.** 1189 / 1213 test spectra are bit-identical to a train row,
  covering every molecule, yet that row's structure is the answer for at most ~46 / 400.
- **Candidate recall.** Whether the answer is among the ~1038 candidates per molecule is unknown;
  rank and recall cannot be separated from one MRR.
- **29 / 400 molecules get no forward-model evidence in our pipeline** (no positive `[M+H]+` / `[M+Na]+`
  spectrum); the forks' coverage of the hidden test is unknown.
- **The image and RDKit switch** (v11 → v12) was not measured on its own; v11's output was lost to
  a session restart before it could be compared.
