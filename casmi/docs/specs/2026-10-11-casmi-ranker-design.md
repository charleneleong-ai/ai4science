# CASMI ranker — design

Replace the hand-set score of the 0.349 notebook (z(FPNet f·z) + 0.25 × popularity, then forward-model
re-ordering) with a LightGBM LambdaRank model trained on features we already compute.

## Why train our own

The public `ranker_v4b.pkl` reads 160 features from the public engine (library, analogue, fragment,
subformula channels, five feature families, three extra networks). Porting it means running that engine,
and its library features point at the copied train structures this work showed are almost never the
answer. A ranker over our own features keeps every input explainable.

## Labels — honest by construction

Training queries are train molecules from the
[fold-safe split](https://www.kaggle.com/datasets/dmitriigluzdov/casmi26-fold-safe-fpnet) (`split.parquet`,
`fold == "hold"`, 55,102 structures), scored with `fpnet_fold0.pt`, which never saw them. The label is
1 for the candidate whose InChIKey14 is the molecule's own, else 0. Using `fpnet_full1` here would leak:
it was trained on these molecules and ranks them near-perfectly (wiring check: 0.943).

At test time the same features come from the production FPNet ensemble. The shift between the fold0
checkpoint (training features) and the ensemble (test features) is a known risk; z-scoring each spectral
score within the molecule limits it.

## Candidates and features

Per query: the production PubChem shortlist — the 1000 most-documented connectivities in 5 ppm. **The
ranker ranks only this shortlist, at training and at test time.** The train window is left out of both:
a train query's own structure always sits in it, so training on it would teach "prefer train
structures", and at test time those are almost never the answer (v11). The copied-structure exclusion
then has nothing left to remove.

| feature | source |
|---|---|
| `fz_z`, `fz_rank` | FPNet f·z, z-scored and ranked within the molecule |
| `lsid`, `lpmid`, `lpatent`, `pop`, `pop_rank` | popularity prior |
| `np_coconut`, `np_lotus`, `np_npatlas` | natural-product flags (`pc_np` bits) |
| `ppm_error` | candidate mass vs neutral precursor mass |
| `n_isomers`, `n_candidates` | shortlist members sharing the formula; shortlist size |

Left out, and why:
- **raw `fz`:** training scores come from `fpnet_fold0`, test scores from the two-checkpoint ensemble,
  so the raw scale shifts; only within-molecule z-scores and ranks transfer.
- **fragment explainability:** about 2M candidates make it too slow for one Kaggle run, and the
  notebook does not compute it once FPNet is attached.
- **`in_train`:** it would encode the label at training time.

`handset` (today's z(fz) + 0.25 × popularity) is stored beside the features as the baseline to beat,
not trained on.

Forward-model scores (ICEBERG, GLACIER) stay a separate re-ordering stage after the ranker: computing
them for training queries needs GPU time we can't spend on thousands of molecules.

## Evaluation and the gate

- **Local, valid for this purpose:** MRR@25 on held-out *query molecules* (grouped split of the hold
  fold), ranker vs the current hand-set score, on identical candidate lists. This measures ranking on
  train-like molecules only; it is a go / no-go for submitting, not an estimate of the leaderboard.
- **Gate G1, unchanged:** the submission must beat 0.349 by more than the ±0.005 noise, with the forward
  models applied after it exactly as in v15.

## Plan

1. Training-set builder: a Kaggle kernel that writes features + labels for a sample of hold-fold
   molecules (PubChem tier and prior are mounted there, not local). Start with 2000 molecules.
2. Local LambdaRank training and grouped evaluation against the hand-set score.
3. If it wins locally, ship the model as a dataset and add a `RANKER` stage to the notebook; one
   submission.

## What this cannot establish

Whether test molecules — whose train-row copies carry the wrong structure — rank like train molecules.
Only the leaderboard measures that.
