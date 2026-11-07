# Seed retrains of the public spectrum model (2026-11-03 to 2026-11-06)

**Result:** blending the public scorer with three seed retrains of it scored **0.463** against 0.479 (fork479 kernel
version 2, submission 58046573, 2026-11-06): DISCARD, −0.016. The same-model ensemble route is closed.

## What was trained

sanpier's public kit `casmi_spectrum_model_v1` (Saner Turhaner's spectrum encoder + fingerprint MLP; width 192,
3 layers, 2048-bit fingerprint), read in full before it ran. [`run_seed.py`](../../scripts/seed_retrain/run_seed.py)
copies the kit, applies the shipped `normalized_smiles` alias fix to `prepare.py`, and changes only `seed`; the
schedule is the kit's own (16 × 500 steps holdout, then a 16-epoch full refit), because longer schedules scored lower
in public (40 × 2000 steps: 0.453 against 0.472; 80 epochs: 0.460).
[`launch.sh`](../../scripts/seed_retrain/launch.sh) starts each seed detached on the shared A100 with an
`autoresearch-heartbeat` beside it.

| seed | holdout diagnostic MRR@25 | wall time |
|---|---|---|
| shipped model | 0.484 | — |
| 20261101 | 0.481 | ~1 h 26 min |
| 20261102 | 0.487 | ~1 h 26 min |
| 20261103 | 0.465 | ~1 h 26 min |

Each seed draws a different holdout, so these numbers are not paired. Data: 2,430,328 spectra, 273,833 molecules
(matching the kit's audit). Prep took ~66 min on CPU, training ~20 s per epoch. The three `final_model.pt` files are
in the private dataset `charyeezy/casmi26-seed-scorers`.

## How they were used

The notebook cell in [`fork479.ipynb`](../../notebooks/fork479/fork479.ipynb) ([`b9627b0`](https://github.com/charleneleong-ai/ai4science/commit/b9627b0))
scores the exported candidates with each seed, then replaces the ranks 2–25 neural score with the mean per-molecule
z-score of the four scorers (equal weights), averages the other neural columns, and keeps a neural score only where
all four scored. Any failure restores the public file and the gate rejects the run.

Run evidence: 400 molecules, 688 rows, 300 with every member scoring; Spearman to the public scorer 0.913, 0.933,
0.928; gate passed.

## Not established

- Why it lost. The seeds rank close to the public scorer and two of three validate lower, so equal weights may dilute
  it; a public-weighted blend (e.g. 0.5 public + 0.5 mean of seeds, as in prvsiyan's unrun public notebook) was not
  submitted, since it is a smaller version of the same change.
- Whether seeds from a different recipe or architecture would help as ensemble members.
