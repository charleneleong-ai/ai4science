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

Per query: the production candidate set — the 1000 most-documented PubChem connectivities in 5 ppm plus
train structures in the window — so the ranker sees the same distribution it will rank at test time.
The copied-structure exclusion is **not** applied to training queries: a train query's spectrum is its
own train row, so the exclusion would remove the answer itself. At test time it still runs before the
ranker, as in v11–v15.

| feature | source |
|---|---|
| `fz`, `fz_z`, `fz_rank` | FPNet f·z, raw / z-scored / rank within the molecule |
| `lsid`, `lpmid`, `lpatent`, `pop`, `pop_rank` | popularity prior |
| `np_coconut`, `np_lotus`, `np_npatlas` | natural-product flags (`pc_np` bits) |
| `ppm_error` | candidate mass vs neutral precursor mass |
| `explained` | fragment explainability (existing code) |
| `n_isomers` | candidates sharing the formula in the window |
| `in_train` | candidate is a train structure |

Forward-model scores (ICEBERG, GLACIER) stay a separate re-ordering stage after the ranker: computing
them for training queries needs GPU time we can't spend on 55k molecules.

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
