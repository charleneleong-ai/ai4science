# Experiment card — CASMI 2026 molecule identification

Pattern: [*A generated protein is not a result*](https://rewirebio.io/blog/a-generated-protein-is-not-a-result/),
applied to a retrieve-and-rank campaign. Status as of 2026-10-11; the code is
[`notebooks/kaggle_submission.py`](../notebooks/kaggle_submission.py), the history
[PR #80](https://github.com/charleneleong-ai/ai4science/pull/80).

## Target — frozen

| | |
|---|---|
| task | rank up to 25 candidate structures (SMILES) per molecule from its LC-MS/MS spectra |
| metric | MRR@25 on the tautomer-canonical InChIKey connectivity block (first 14 characters) |
| test set | `test.parquet`: 400 molecules, 1213 spectra; answers hidden |
| scored on | the **public** leaderboard split only; the private split is not visible until the competition closes |
| rules | code competition: notebook only, internet off; about 5 submissions a day |

## Inputs

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

- **19 submissions, 19 scored, 0 rejected** (one an identical re-run of v19 to measure noise: 0.353 vs 0.356). Every one is in the
  [score table](https://github.com/charleneleong-ai/ai4science/pull/80), including the four not
  kept (0.049, 0.254, 0.345, 0.346), plus v18 (0.353, within noise) v19 (0.356, kept), v20 (0.354) and v21 (0.333) in
  [PR #85](https://github.com/charleneleong-ai/ai4science/pull/85).
- **Two measurement kernels, not submissions:** the ranker training set (2000 held-out train
  molecules, 3 h 37 min) and the second-pass recall probe (300 molecules, 35 min; v1 replaced before
  it finished because an empty mass window would have crashed it).
- **18 notebook versions run on Kaggle;** the 3 not submitted were setup or reproducibility runs
  (a missing RDKit wheel, a push-from-repo check).
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
| second pass skips candidates with Tanimoto ≥ 0.5 to the copied structure | **after** v18 scored | v18's changed top-1 picks drifted toward the copied structures; only 0.5 tried |

Choices marked **after** were selected on the public leaderboard, so the reported 0.356 is an
optimistic estimate of the private score. Differences under about ±0.005 — the run-to-run noise
of identical runs — were treated as no change.

## What is not measured

- **The private score.** Every number here is public-split.
- **Any local estimate.** Train-derived labels were wrong (a local 0.664 scored 0.115), and
  consensus labels from 10 public notebooks failed calibration: 3 lineages agreeing on 12 / 400
  top-1 picks.
- **Why the data is shaped as it is.** 1189 / 1213 test spectra are bit-identical to a train row,
  covering every molecule, yet that row's structure is the answer for at most ~46 / 400.
- **Candidate recall.** Whether the answer is among the ~1038 candidates per molecule is unknown;
  rank and recall cannot be separated from one MRR.
- **29 / 400 molecules get no forward-model evidence** (no positive `[M+H]+` / `[M+Na]+` spectrum).
- **The image and RDKit switch** (v11 → v12) was not measured on its own; v11's output was lost to
  a session restart before it could be compared.
