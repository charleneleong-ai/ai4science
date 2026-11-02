# CASMI26 leaders sweep, 2026-11-01

This is a delta on [the 2026-10-30 report](../../reports/2026-10-30-casmi-sota-delta.md). It records only what is
new.

Labels:
- **CONFIRMED**: seen directly.
- **HINTED**: someone stated it, but not a leader.
- **SPECULATION**: inference.

Notebooks read (never executed) are in `~/casmi-data/ref/leaders_20261101/`.

## Leaderboard (Kaggle CLI, 2026-11-01)

| Team | Score |
|---|---|
| Savko | 0.521 |
| Oleksandr Parfionov | 0.519 |
| democatXamer (new to the top 5) | 0.512 |
| pikachu | 0.504 |
| Shehab Alshehabi | 0.502 |
| CHU1PC | 0.495 |
| team051, Paweł Biernat | 0.490 |
| JJ, "Day: chem." | 0.489 |

- The top moved from 0.504 (10-25) to 0.517 (10-30) to 0.521 (11-01).
- On 10-25 the leaders had 65 (Parfionov), 167 (Savko) and 128 (pikachu) submissions, according to a
  participant's snapshot ([cosmorin](https://www.kaggle.com/code/cosmorin/why-1-is-1-enveda-casmi26-gap-analysis)).
  These counts are participant-reported.

## Evidence about leaders

- **CONFIRMED: no leader publishes anything.**
  - Accounts checked: parfyon, savchenko, w1326627677 (democatXamer), shehabalshehabi, tadashikato and
    lowasitgoes (CHU1PC; teammate Christoffer Thimsen is a Master), pawebiernatingenix (created 2026-10-12) and
    marvintmb.
  - None has a public dataset or a CASMI notebook.
  - Their CASMI forum activity is Q&A only.
- **CONFIRMED: no public notebook sits between 0.479 and the leaders.** The top of the Code tab is all forks of the
  0.465 → 0.477 → 0.479 lineage.
- **HINTED (non-leader, rank ~2650):** Stefan Dinu declares a harmonised GNPS + MassSpecGym library
  ([754574](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754574);
  [Zenodo 20042904](https://zenodo.org/records/20042904), CC BY 4.0).
  - The record holds 444k spectra and 45k SMILES; he uses 185k spectra covering 38k structures, plus COCONUT.
  - This is a much larger reference library than the 751 external spectra in the lineage.
- **SPECULATION:** cosmorin argues that the gap lies upstream, in candidate recall and the neural model, and that
  the leaders' submission volume points to a private validation harness. This is inference, not evidence.
- **Others have tried transduction:**
  [youyouwou v113transd](https://www.kaggle.com/code/youyouwou/casmi26-v113transd) re-ranks by test-neighbour
  votes (τ 0.70, λ 1.5) and reports no score.

## Host rulings since 10-30

- [753662](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753662)
  (10-30, David Healey): pre-built datasets are allowed "as long as they are free and openly accessible and they
  were not trained on or built from proprietary data".
  - This partly defuses the asset-eligibility risk.
  - Still open: datasets that contain weights (753900 follow-up), test-time weight updates (748882), BMDMS-NP and
    [754610](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754610).
- [745134](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/745134)
  (10-28): the [ICEBERG PubChem atlas](https://iceberg-ms.mit.edu/) may be shipped as a private dataset with its
  retrieval code.
  - The atlas holds predicted spectra for 91.4M PubChem structures, from a NIST'23-trained ICEBERG.
  - A membership-only list is also allowed.
  - **No public notebook uses it**: a `grep` over the 0.479, 0.477 and 7 other current notebooks found 0 hits.
    Our fork does not use it either.

## Candidate experiments (one change each, ranked)

| # | Change | Zero-submission reach check | Gain | Cost |
|---|---|---|---|---|
| 1 | Use the ICEBERG atlas as an in-silico library, scoring PubChem-channel and same-formula candidates by cosine to the atlas spectrum | Measure the atlas's InChIKey14 coverage of Enveda-180/panel golds that sit in the pool below rank 1 (the atlas **excludes NIST'23 structures**, which may remove many documented natural products). Then measure the paired MRR delta on the panel. Submit only if about 8 public-equivalent molecules move and the delta exceeds 2×SE | +0.000-0.015 (speculative) | 1-2 days, plus the atlas download by formula |
| 2 | Extend the reference library with Zenodo 20042904 | Count its InChIKey14s that are absent from the train.parquet structures and from the 751-spectrum table. Count how many Enveda-180 golds become library-matchable at cosine ≥ 0.7 | +0.000-0.010 | 0.5 day |

## Papers

Forum hints point at no particular method. An HF Papers search turned up one item not in the 10-30 report:
[dIon, arXiv 2610.06282](https://arxiv.org/abs/2610.06282) (10-05), self-supervised learning for MS/MS. Nothing
links it to a leader.

## Not found

- Any leader's method.
- Team submission kernels, because the MCP leaderboard and submission calls returned Unauthenticated.
- Per-notebook public scores, which the CLI omits.
- Usernames for pikachu, JJ and team051.
- The Discord.
