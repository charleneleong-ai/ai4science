# Experiment card — CASMI 2026 molecule identification

Pattern: [*A generated protein is not a result*](https://rewirebio.io/blog/a-generated-protein-is-not-a-result/),
applied to a retrieve-and-rank campaign. Status as of 2026-11-13 (AEDT): the current best (0.487, inside the noise of
0.479) is a fork of Hongbin Gu's public 0.487 notebook, unchanged except for a smoke commit run, [`notebooks/fork487/fork487.ipynb`](../notebooks/fork487/fork487.ipynb);
the base before it (0.479, inside the noise of 0.477) is an unchanged fork of alex chilton's public 0.479 notebook, [`notebooks/fork479/fork479.ipynb`](../notebooks/fork479/fork479.ipynb);
the base before it (0.477) is an unchanged fork of Busya PRIME's public 0.477 notebook, [`notebooks/fork477/merged_fork.ipynb`](../notebooks/fork477/merged_fork.ipynb);
before it (0.465) is an unchanged fork of dmitriigluzdov's public notebook at version 26, [`notebooks/fork464/data_fork.ipynb`](../notebooks/fork464/data_fork.ipynb);
before it (0.456), an unchanged fork of Hongbin Gu's public 0.456 notebook, [`notebooks/fork456/pubchemlite_fork.ipynb`](../notebooks/fork456/pubchemlite_fork.ipynb);
before it (0.446), an unchanged fork of the public E006 notebook,
[`notebooks/fork446/e006_fork.ipynb`](../notebooks/fork446/e006_fork.ipynb); the first fork (0.415) is of the 0.421 notebook,
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

## Current best — extlib-fork, public MRR@25 0.487 (PROMOTE, inside the noise)

| | |
|---|---|
| kernel | `charyeezy/casmi-fork487` (private), kernel version 1, submission 58272211 (2026-11-12 23:06 AEDT, scored by 03:43 AEDT 13 November), commit [`2dc7d19`](https://github.com/charleneleong-ai/ai4science/commit/2dc7d19); reproduces the source's public 0.487; +0.008 over dedup-fork, one to two molecules of the ~124-molecule public split and inside the ±0.016 between submissions, so PROMOTE as the new base, not a measured gain |
| source | [casmi26-v589-alex0484-extlib](https://www.kaggle.com/code/guhongbin/casmi26-v589-alex0484-extlib) (Hongbin Gu), run unchanged by huichieh for 0.487 (submission 58185447, 2026-11-10); md5 of the source `ec45cc46672d00245bdb6bb0c39712f7`; one change, `COMMIT_MODE = 'smoke'` (was `'fast'`), so the commit run executes every stage on the 12 visible molecules; the scoring rerun is full either way. See [`notebooks/fork487/README.md`](../notebooks/fork487/README.md) |
| method | popwenman's Gold Apex V22, with: the public spectrum model and an 80-epoch retrain (popwenman's public `casmi26-scorer-train-v80` kernel) combined as 0.5·z(public) + 0.5·z(ours); romanrozen's timsTOF mass prior; a documentation prior on ranks 2–25 with a LOTUS backfill; external-library rank-1 promotion (cosine ≥ 0.8, K = 4) against Takumu Hata's extlib (200,724 spectra) plus harmonised MassBank, 324,884 in all |
| datasets | 10 plus one kernel source, as listed in [`kernel-metadata.json`](../notebooks/fork487/kernel-metadata.json); the source's `id_no` and `keywords` removed (a copied `id_no` overrides `id` and makes the push fail with 403) |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461`, NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_fork487.py` (outside the repo) passed on the commit run: the fork477 gate, plus `OURS BUNDLE READY`, the ensemble line, mass prior on (394 / 400 lists touched, 447 candidates penalised), documentation prior reordered 12 rows, extlib loaded 324,884 spectra with 0 matches and 0 promotions (the smoke molecules are train copies), no fallback or unexpected `Note: ... exception` line. Outputs: `~/casmi-data/kernel_outputs/fork487-v1` |
| next | kernel version 2 ([`28c5ce3`](https://github.com/charleneleong-ai/ai4science/commit/28c5ce3)), the learned isomer reranker on (`V8 = 'BOTH'`, `RR8_ALPHA = 0.20`), gate `gate_fork487_v2.py`: 0.487 (58292050), no change, discarded; kernel version 3 ([`62b95fb`](https://github.com/charleneleong-ai/ai4science/commit/62b95fb)), `RR8_ALPHA = 0.15`, running |

The settings were fixed by the notebook's authors before we ran it; we chose none of them.

## Kept — dedup-fork, public MRR@25 0.479 (KEEP, inside the noise)

| | |
|---|---|
| kernel | `charyeezy/casmi-fork479` (private), kernel version 1, submission 57867559 (2026-11-01 20:37 AEDT), commit `a719425`; reproduces the source's public 0.479; +0.002 over merged-fork, below G1's +0.005 and inside the source's own run-to-run spread (~0.002), so KEEP: adopted as the base, not a measured gain |
| source | [CASMI26 LB 0.479 \| MIDX Exact Dedup Fill](https://www.kaggle.com/code/alexchilton/casmi26-lb-0-479-midx-exact-dedup-fill) (alex chilton), which runs Baidalin Adilzhan's *casmi26 midx exact dedup fill* unchanged; byte-identical to its version as of 2026-11-01 (md5 `a6871ca326db69f8c4a5caff8d6b7d36`); only [`kernel-metadata.json`](../notebooks/fork479/kernel-metadata.json) differs (our id, private). Vendored in [`notebooks/fork479/README.md`](../notebooks/fork479/README.md) |
| method | the merged-fork pipeline below, then a final layer on ranks 2–25: a documentation prior from `pool_popularity.csv` (PubChem popularity 0.10, spectral-library bonus 0.05 only where the best library match reaches 0.87, natural-product bonus 0.02); the old leader pinned at rank 2 after a fragment-check move; deduplication by the metric's InChIKey14 key, freed slots filled with scored outside candidates. Rank 1 is left to the merged-fork logic |
| datasets | 8, as listed in [`kernel-metadata.json`](../notebooks/fork479/kernel-metadata.json): merged-fork's seven plus `metric/rdkit-2026-3-3-wheel` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461`, NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_fork479.py` (outside the repo) passed on the commit run: the merged-fork gate, plus the final layer's log lines: the prior loaded 710,701 pool records and reordered ranks 2–25 in 12 / 12 smoke rows, 0 dedup swaps and 0 outsider fills on the smoke run, no `Note: ... exception` line (no neutral-prior fallback taken). Outputs: `~/casmi-data/kernel_outputs/fork479-v1` |

The settings were fixed by the notebook's authors before we ran it; we chose none of them. The source's own table
shows the same code about 0.002 apart between runs, and two notebooks reported at 0.478 rescored 0.468 and 0.476, so
0.479 is numerically our highest public score but not distinguishable from 0.477.

## Discarded — mass-prior, public MRR@25 0.479 (fork479 kernel version 3)

| | |
|---|---|
| change | romanrozen's timsTOF precursor-mass prior from popwenman's public [Gold Apex V22](https://www.kaggle.com/code/popwenman/casmi26-0482-gold-apex-v22), alone on the v1 notebook, constants unchanged ([`5f8cfb6`](https://github.com/charleneleong-ai/ai4science/commit/5f8cfb6)): ppm = (exact − neutral) / neutral, neutral = median precursor m/z − adduct; ranks 2–25 lose 0.15 · min(1, (\|ppm + 1.5\| − 3) / 5) beyond 3 ppm; an in-band rank 2–3 takes rank 1 from a leader more than 5 ppm off when the best library match is below 0.9. V22's scorer blend, prior weights and LOTUS backfill left out |
| checked before | on 2,138 timsTOF train spectra the true structure sits at median −1.9 ppm, 99.3% within 3 ppm of −1.5, 0.14% beyond 5; on the 12 smoke lists of v1's commit run 59 of 288 rank 2–25 candidates fall outside the band (6 / 12 lists), the rank-1 rule moved none |
| run | submission 58107187 (2026-11-08 10:29 AEDT), scored by 14:52 AEDT; gate `gate_fork479_v3.py` passed (fork479 gate, `MASS PRIOR` on with lists touched, `MASS LEAD` ran, no fallback); rank-1 rule 0 moves on the commit run |
| result | 0.479, no change against 0.479: DISCARD |
| fixed before seeing the score | V22's constants, the gate, the bar (≥ 0.485) |

## Discarded — seed-blend, public MRR@25 0.463 (fork479 kernel version 2)

| | |
|---|---|
| change | on ranks 2–25, the neural score is the mean per-molecule z-score of four scorers: the public v1 spectrum model and three retrains of it; the other neural columns are averaged and a candidate keeps a neural score only where all four scored it ([`b9627b0`](https://github.com/charleneleong-ai/ai4science/commit/b9627b0)) |
| retrains | sanpier's public kit `casmi_spectrum_model_v1` read in full, run unchanged except the shipped `normalized_smiles` alias fix and the seed (20261101–3); 16 × 500 steps holdout then a 16-epoch refit, on the shared A100 (about 1 h 26 min each, 2,430,328 spectra, 273,833 molecules); holdout 0.481 / 0.487 / 0.465 against the public model's 0.484; private dataset `charyeezy/casmi26-seed-scorers` |
| run | submission 58046573 (2026-11-06 16:14 AEDT), scored by 21:03 AEDT; gate passed (all three seeds ran, one `SEED BLEND applied` line, no fallback); 400 molecules, 688 exported rows, 300 with every member scoring; Spearman to the public scorer 0.913 / 0.933 / 0.928 |
| result | 0.463, −0.016 against 0.479, beyond the ±0.010 noise: DISCARD |
| fixed before seeing the score | equal weights, the shipped schedule (longer schedules had scored 0.453 and 0.460 in public), the gate |

## Previous best — merged-fork, public MRR@25 0.477

| | |
|---|---|
| kernel | `charyeezy/casmi-merged-fork` (private), kernel version 1, submission 57822991 (2026-10-31), commit `c467db3`; reproduces the source's public 0.477; +0.012 over data-fork, above the ±0.005 noise (G1: PROMOTE) |
| source | [LB 0.477 · Two 0.472 Fixes, Merged the Right Way](https://www.kaggle.com/code/busyaprime/lb-0-477-two-0-472-fixes-merged-the-right-way) (Akmal Xodarev, Busya PRIME), byte-identical to its current version as of 2026-10-28 (md5 `74efef20daa877c43243024d830ef032`), the version its page reports at public 0.477; only [`kernel-metadata.json`](../notebooks/fork477/kernel-metadata.json) differs (our id, private). Vendored with its version table in [`notebooks/fork477/README.md`](../notebooks/fork477/README.md) |
| method | runs the pinned 0.465 data-fork pipeline once, then stacks two public changes that each showed 0.472: (1) alex chilton's leader-safe neural re-rank, reciprocal-rank fusion of each molecule's final 25 with the spectrum model at weight 0.25 for rank 1 and 0.75 for ranks 2–25; (2) Roman Tamrazov's Spectral Parallax fragment check, which may replace rank 1, its leader moves taken from the 0.50 (V2) order and merged as `[P0, V0]` + the rest of the CUSTOM order. Candidate sets never change |
| datasets | 7, as listed in [`kernel-metadata.json`](../notebooks/fork477/kernel-metadata.json): `sanpier/casmi-spectrum-model-v1`, `sanpier/catalogue-experiment-saner`, `dmitriigluzdov/casmi26-{fold-safe-fpnet, natural-product-knowledge-table, natural-product-panel-benchmark, pubchem-popularity-prior, attributed-v4g-reproduction-assets}` (the same seven as neural-fork) |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:9a6bcfda6d03b2cd9aa347ee8c1bae62043328ed0f4afe9bda014e4f8bbf9435` (the source's own GPU image, newer than the one our other forks pin), NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_fork477.py` (outside the repo) passed on the commit run: the 0.465 pipeline's checks from its stage logs, both layer reports written, the fragment check not skipped, 400 rows with the 12 smoke molecules full |

The settings were fixed by the notebook's authors before we ran it; we chose none of them. Both stacked
changes were scored on the public leaderboard upstream (0.472 each, 0.477 merged), so 0.477 may be
optimistic on the private split.

**Discarded — analogue-fork, public MRR@25 0.477** (fork477 v3, kernel version 2, submission 57826104, commit `cfabe09`). One cell after the
merge adds analogue propagation across test molecules, allowed by the host on 2026-10-26 (discussion
752297): molecules whose leader has a library match ≥ 0.87 are anchors; for every other molecule whose
spectra reach cosine ≥ 0.70 with an anchor, candidates below rank 1 gain 0.10 × cosine × Tanimoto(candidate,
anchor leader). Rank 1 and the candidate set never change. At most 113 / 400 visible-test molecules have
such a neighbour. Fixed **before** submission. Gate `gate_fork477_v3.py` passed: the stage applied, and no
leader or candidate set changed against `submission_pre_analogue.csv`. It scored 0.477 (2026-10-31),
identical to the base: no measurable effect at λ 0.10, discarded. The hidden rerun's logs are not visible,
so whether few lists moved or the moves did not matter is not established.

## Earlier best — data-fork, public MRR@25 0.465

| | |
|---|---|
| kernel | `charyeezy/casmi-data-fork` (private), kernel version 1, submission 57588653 (2026-10-24); +0.009 over pubchemlite-fork, above the ±0.005 noise (G1: PROMOTE) |
| source | [CASMI 26: The Data Behind Spectra to Structures](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-the-data-behind-spectra-to-structures) (dmitriigluzdov), version 26, byte-identical as of 2026-10-24 (md5 `4edf7f5dbb3bee2aff0cc52f63961388`); only [`kernel-metadata.json`](../notebooks/fork464/kernel-metadata.json) differs (our id, private). Vendored with its version table in [`notebooks/fork464/README.md`](../notebooks/fork464/README.md) |
| lineage | the 0.456 PubChemLite recipe (our pubchemlite-fork) → GLACIER proxies kept where ICEBERG scores exist, joint ICEBERG+GLACIER re-ordering of proxy-covered lists (0.464) → five fold-safe FPNets in the PubChem channel plus a sixth PubChem slot at rank 24 (0.464) → candidate-aware rank-1 PubChem promotion (Saner Turhaner, 0.465 upstream). Kaggle serves only version 26, which was not scored upstream; we reproduced only the endpoint |
| datasets | 5 of the author's datasets attached unchanged, as listed in [`kernel-metadata.json`](../notebooks/fork464/kernel-metadata.json): `dmitriigluzdov/casmi26-{fold-safe-fpnet, natural-product-knowledge-table, natural-product-panel-benchmark, pubchem-popularity-prior, attributed-v4g-reproduction-assets}` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461` (the same GPU image as the other forks), NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_fork451.py` (outside the repo) passed on the commit run: output 400 rows; 12 smoke molecules full; ICEBERG and GLACIER scored 11 each |

The settings were fixed by the notebook's authors before we ran it; we chose none of them. Its lineage
was scored on the public leaderboard upstream, so 0.465 may be optimistic on the private split.

## Sixth fork — neural-fork, public MRR@25 0.444 (discarded)

| | |
|---|---|
| kernel | `charyeezy/casmi-neural-fork` (private), kernel version 1, submission 57615727 (2026-10-25); [`notebooks/fork466/neural_fork.ipynb`](../notebooks/fork466/neural_fork.ipynb), vendored in commit `2a2411c`, edited in `9ae5247` |
| source | [PubChem + Neural025 — LB 0.466](https://www.kaggle.com/code/sanpier/pubchem-neural025-lb-0-466) (Saner Turhaner), vendored byte-identical to its current version as of 2026-10-25 (md5 `af68cfeeb06a251186cf043684ca3be2`, before our edit); [`kernel-metadata.json`](../notebooks/fork466/kernel-metadata.json) sets our id, private, and drops an empty dataset entry. Vendored with its version table in [`notebooks/fork466/README.md`](../notebooks/fork466/README.md) |
| method | runs a frozen 0.465 baseline, then re-orders each molecule's final 25 by reciprocal-rank fusion `1/(3 + baseline_rank) + w/(3 + neural_rank)` with a separately trained spectrum model (`sanpier/casmi-spectrum-model-v1`); candidate sets never change. Upstream `w = 0.25` scored 0.466; the current version carries 0.15, not scored upstream |
| change | one: `NEURAL_WEIGHT` 0.50 (cells 1 and 9), fixed **before** submission |
| datasets | 7, as listed in [`kernel-metadata.json`](../notebooks/fork466/kernel-metadata.json): `sanpier/casmi-spectrum-model-v1`, `sanpier/catalogue-experiment-saner`, `dmitriigluzdov/casmi26-{fold-safe-fpnet, natural-product-knowledge-table, natural-product-panel-benchmark, pubchem-popularity-prior, attributed-v4g-reproduction-assets}` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7d…a461` (the same GPU image as the other forks), NVIDIA T4, internet off |
| run evidence | gate `~/casmi-data/gates/gate_neural.py` (outside the repo) passed on the commit run: output 400 rows; 12 smoke molecules full; ICEBERG and GLACIER scored 11 each; the neural stage at weight 0.5 re-ordered all 12 smoke molecules and changed top-1 on 2 |
| result | 0.444: −0.022 against `w = 0.25` upstream and −0.021 against data-fork's 0.465, which stays the best; discarded. The neural re-rank adds at most +0.001 to the 0.465 base; more weight hurts |

## Earlier best — pubchemlite-fork, public MRR@25 0.456

| | |
|---|---|
| kernel | `charyeezy/casmi-pubchemlite-fork` (private), kernel version 1, submission 57435777 (2026-10-20); reproduces the source's public 0.456 exactly; +0.010 over e006-fork |
| source | [casmi26-v465-v15-pop015](https://www.kaggle.com/code/guhongbin/casmi26-v465-v15-pop015) (Hongbin Gu), byte-identical as of 2026-10-20 (md5 `9d1ab164357187405aa67fd1335d081d`); only [`kernel-metadata.json`](../notebooks/fork456/kernel-metadata.json) differs (our id, private). Vendored with its lineage table in [`notebooks/fork456/README.md`](../notebooks/fork456/README.md) |
| lineage | the 0.451 base (E006 + `rel > 150` + `rank2` hedge + standardised ChEBI/LIPID MAPS candidates) → PubChemLite structures in engine 2's library (0.455) → N1A: GLACIER proxy scores from other adducts for molecules without an `[M+H]+` spectrum → PubChem popularity weight 0.25 → 0.15 (0.456). Step scores are the author's public ones; we reproduced only the endpoint |
| ranker | as e006-fork: engine 1 `ranker_v4b.pkl` (160 features, 4 boosters); engine 2 rankers fitted in the notebook; none retrained by us |
| datasets | 15 public datasets attached unchanged, as listed in [`kernel-metadata.json`](../notebooks/fork456/kernel-metadata.json), including `franciscoangulo/casmi26-pubchemlite-cc0-candidate-fingerprints` (the PubChemLite candidates), `sergiomalv/casmi-chebi-lmsd-structures`, `prvsiyan/coconut-casmi26-candidates`, `megayak/casmi26-simulated-ranker-rows`, `ahmedberatozer/casmi26-{fpnet-full1, glacier, iceberg, pubchem-tier, v2-pool, v3-models, v4b-models}` and `metric/rdkit-2026-3-3-wheel` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7dd9c54116ecd1bcc88817c5469b88387388fade02bfa8bf3fc647d461` (the same GPU image as the other forks), NVIDIA T4, internet off |
| run evidence | gate passed on the commit run (about 1190 s): ranker 160 features / 4 boosters; errors 0; ICEBERG and GLACIER scored 11 smoke molecules each; output 400 rows |

The settings were fixed by the notebook's authors before we ran it; we chose none of them. The
popularity weight 0.15 was chosen on the public leaderboard upstream, so 0.456 may be optimistic on the
private split.

## Earlier best — e006-fork, public MRR@25 0.446

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
optimistic on the private split. One setting has been changed since: kernel version 3 (e006-v3)
set the engine-2 weight in the final rank fusion, `ALPHA` in cell 35, to 0.8 instead of 0.6 (a one-byte
change); its gate passed (ranker 160 / 4, errors 0, ICEBERG and GLACIER 11 smoke molecules each) and it
scored 0.437 (submission 57310111), −0.009: discarded. The
[2026-10-17 SOTA sweep](../reports/2026-10-17-casmi-sota-sweep.md) then named the host-approved FRIGID de novo
model as the next lever; it was closed offline (below).

## Third fork — spectra-fork, public MRR@25 0.445 (discarded)

| | |
|---|---|
| kernel | `charyeezy/casmi-spectra-fork` (private), submission 57311145; [`notebooks/fork451/spectra_fork.ipynb`](../notebooks/fork451/spectra_fork.ipynb) |
| source | [CASMI 26: From Spectra to Structures](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-from-spectra-to-structures) (dmitriigluzdov), forked unchanged at its current version 53 (md5 `951a34f88affa3b89d7fdefc166e60f5`). Its own table reports 0.451 at version 52 (E006 + depth-7 engine-2 rankers + standardised ChEBI/LIPID MAPS candidates); Kaggle serves only version 53, which adds ranker seeds |
| datasets | the author's 4, as listed in [`kernel-metadata.json`](../notebooks/fork451/kernel-metadata.json): `dmitriigluzdov/casmi26-{attributed-v4g-reproduction-assets, natural-product-knowledge-table, natural-product-panel-benchmark, pubchem-popularity-prior}` |
| runtime | image `gcr.io/kaggle-private-byod/python@sha256:37c64f7d…a461` (GPU), NVIDIA T4, internet off |
| result | 0.445, inside the noise of 0.446 and 0.006 below the reported 0.451; discarded. Whether version 53's seeds or non-reproduction explain the gap is not established |

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

## Offline experiments — no submission

**MIST `msg_aug` encoder probe, FAIL (2026-11-01).** Run on an A100 80 GB, about 21:30 AEDT; notes:
[`research_notes/2026-11-01-mist-msg-aug/notes.md`](../research_notes/2026-11-01-mist-msg-aug/notes.md).

| | |
|---|---|
| models | `mist_msg_aug.ckpt` and `mist_msg_noaug.ckpt` from [Zenodo 23303289](https://zenodo.org/records/23303289) (md5 `3e921792…` and `80594e4d…`, both matching the record); [coleygroup/MIST](https://github.com/coleygroup/MIST) at `48e339be` |
| sample | the same 50 Enveda-180 FRIGID targets as below, plus the two 40-spectrum check samples; 2 of 130 spectra had no annotated peaks and stay in the denominator |
| gate | mean fingerprint Tanimoto ≥ 0.40 on the 50 FRIGID molecules (4096-bit Morgan r2), **fixed before the run** |
| result | `msg_aug` 0.385 (median 0.336, SE 0.022, 17 / 50 at ≥ 0.4), against 0.261 for the old `mist_msg.pt`: 0.015 short, so FAIL; the MIST/FRIGID line stays closed. Reopening it would be a post-hoc decision |
| control | `noaug` is unusable: 0.19 even on its own training spectra (`aug`: 0.86), saved at step 30,000 against 430,965; how much of the gain is the synthetic data rather than code or training length is not measured |

**External GNPS + MassSpecGym library, reach 0 (2026-11-01).** [Zenodo 20042904](https://zenodo.org/records/20042904)
(CC BY 4.0) holds 41,210 structures; 3,096 of them (13,814 spectra) are not already in `train.parquet`, which holds
221k GNPS, 102k MassBank and 92k MoNA spectra. Those 3,096 match **0 / 400** test molecules at cosine ≥ 0.5 within
20 ppm in the same ion mode, so it was dropped without a submission. Caveat: every visible test spectrum matches train
at ≥ 0.87, so the visible test looks like train stand-ins rather than the hidden test; for novel Enveda-180-like
structures exact reach is about 0 by construction.

The two below ran on 2026-10-18 on an Apple M5 (PyTorch 2.6.0, MPS), outside Kaggle; nothing was submitted.
Full record: `~/casmi-data/frigid/RESULTS.md` (outside the repo).

**FRIGID de novo, NO-GO.**

| | |
|---|---|
| models | FRIGID-base: `mist_msg.pt` encoder (MIST-MSG) + `DLM.ckpt` decoder from the official pretrained checkpoints, config `spec2mol_benchmark_msg.yaml` (FP threshold 0.187, up to 10 formula-matched candidates in 100 attempts); a substitute NGBoost token-length model we trained (the release ships none), with the 50 targets removed |
| sample | 50 Enveda-180 molecules with an `[M+H]+` spectrum and no spectra in any other library, drawn with `numpy.default_rng(0)` from 165,122 eligible; 0 excluded; the most-peaked spectrum of each |
| bar | top-10 exact (InChIKey-14) ≥ 5 / 50 with the **true** formula, **fixed before the run** |
| result | 0 / 50 at top 1, 5 and 10; mean max Tanimoto of the top 10 to the answer 0.230; 50 / 50 spectra returned 10 candidates |
| validity | 3588 / 3819 generations valid (94%, a lower bound: early stopping counts unchecked samples of the last batch as invalid); 35 s per spectrum |
| diagnostics, **post hoc** | MIST predicted-fingerprint Tanimoto 0.26 on the 50 (0.25 on 40 random Enveda-180 spectra) against 0.63 on 40 MassBank spectra; the decoder fed the true fingerprint recovers 49 / 50 at top 1, inflated because its 230M-structure pretraining probably includes these compounds. The encoder is the bottleneck |
| declared deviations | the substitute token model; subformulae assigned at 20 ppm (repo default 15) for timsTOF; MIST-MSG chosen over MIST-CANOPUS before the run; the predicted-formula step skipped once the true-formula bound scored 0 |

**FPNet as FRIGID's encoder, closed.** The fold-safe FPNet (`dmitriigluzdov/casmi26-fold-safe-fpnet`,
`fpnet_fold0.pt`) was checked as a replacement for MIST. Only 10 of the first 50 molecules are held out
of its training, so a new paired sample was drawn: 50 of the 32,993 eligible molecules in its held-out
fold, none overlapping the first 50, with MIST re-run on the same spectra. The gates for running the
FRIGID decode were set in the check's protocol: FRIGID's 4096-bit radius-2 Morgan space and mean
Tanimoto ≥ 0.45. Both failed:

| fingerprint vs truth, n = 50 held out | mean Tanimoto |
|---|---:|
| MIST-MSG, 4096-bit Morgan (FRIGID's space) | 0.268 |
| FPNet, own 10,226-bit space | 0.245 |
| FPNet mapped to 4096-bit Morgan (only 1,464 of 4,096 ECFP4 bits kept) | 0.180 |

FPNet beat MIST on 4 / 50 molecules in the shared space. **Post hoc**, after those numbers: the best
global thresholds (0.254 own space, 0.208 mapped, MIST 0.269), the oracle top-k bits, the trained-on vs
held-out comparison (0.235 vs 0.245), and a 50-way retrieval in which FPNet ranks the truth first for 42 / 50
against MIST's 12 / 50, showing it learned to rank rather than to reconstruct bits. Retraining a
4096-bit head was estimated at 10–15 T4 hours and not run; the decode step was never reached.

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
| spectrum model (neural re-rank, merged-fork) | `sanpier/casmi-spectrum-model-v1`, a separately trained spectrum model; RRF weight 0.25 for rank 1, 0.75 for ranks 2–25 | [`sanpier/casmi-spectrum-model-v1`](https://www.kaggle.com/datasets/sanpier/casmi-spectrum-model-v1), attached at its latest version |
| fragment check (merged-fork) | Spectral Parallax (Roman Tamrazov): untrained, MetFrag-style fragment re-scorer, code inside the notebook; may replace rank 1 | [`notebooks/fork477/merged_fork.ipynb`](../notebooks/fork477/merged_fork.ipynb) |
| RDKit | 2026.3.3 in the notebook; 2025.03.6 inside the forward-model subprocesses (their requirement) | [`metric/rdkit-2026-3-3-wheel`](https://www.kaggle.com/datasets/metric/rdkit-2026-3-3-wheel); ICEBERG's wheels |
| runtime | Python 3.12, image `gcr.io/kaggle-private-byod/python@sha256:37c64f7d…a461`, NVIDIA T4 | [`kernel-metadata.json`](../notebooks/kernel-metadata.json) |

FPNet full1 was **trained on all of train**, including the train rows the test spectra copy (see
below). No model here was trained or fine-tuned by this work; every weight is a public checkpoint, and the
fragment check has none.

## Budget — every attempt in the denominator

- **34 submissions, 34 scored, 0 rejected.** 19 from our own pipeline, one an identical re-run of v19
  to measure noise (0.353 vs 0.356); every one is in the
  [score table](https://github.com/charleneleong-ai/ai4science/pull/80), including those not kept
  (0.049, 0.254, 0.345, 0.346), plus v18 (0.353, within noise), v19 (0.356, kept), v20 (0.354) and v21 (0.333) in
  [PR #85](https://github.com/charleneleong-ai/ai4science/pull/85). Two from the 0.421 fork: fork-v1
  (0.415, 57288940) and fork-v2 (0.415, 57290631, discarded). Two from the E006 fork: e006-fork
  (0.446, 57298936, kernel version 2) and e006-v3 (0.437, 57310111, kernel version 3, discarded). One
  from the 0.451 notebook: spectra-fork (0.445, 57311145, discarded). One from the 0.456 notebook:
  pubchemlite-fork (0.456, 57435777, kernel version 1, promoted). One from the version-26 notebook:
  data-fork (0.465, 57588653, kernel version 1, promoted). One from the 0.466 neural re-rank notebook:
  neural-fork (0.444, 57615727, kernel version 1, discarded). One from the 0.477 notebook: merged-fork
  (0.477, 57822991, kernel version 1, promoted); its kernel version 2, analogue-fork (0.477, 57826104, discarded).
  Three from the 0.479 notebook: dedup-fork (0.479, 57867559, kernel version 1, kept within the noise), its kernel
  version 2, seed-blend (0.463, 58046573, discarded), and kernel version 3, mass-prior (0.479, 58107187, discarded). One from the 0.487 notebook: extlib-fork (0.487,
  58272211, kernel version 1, promoted) and isomer-rerank (0.487, 58292050, kernel version 2, discarded).
- **Two measurement kernels, not submissions:** the ranker training set (2000 held-out train
  molecules, 3 h 37 min) and the second-pass recall probe (300 molecules, 35 min; v1 replaced before
  it finished because an empty mass window would have crashed it).
- **18 notebook versions of our own pipeline run on Kaggle,** plus 2 of the 0.421 fork, 3 of the E006 fork, 1 of the 0.451 notebook, 1 of the 0.456 notebook, 1 of the version-26 notebook, 1 of the 0.466 neural re-rank notebook, 2 of the 0.477 notebook and 2 of the 0.479 notebook; of our own pipeline's, the 3 not submitted were setup or reproducibility runs
  (a missing RDKit wheel, a push-from-repo check). The E006 fork's version 1 ran on the CPU image
  its source's metadata names, failed the notebook's GPU preflight and was not submitted.
- **Three spectrum-model retrains** on the shared A100 (seeds 20261101–3, about 4 h 20 min in all, plus a smoke
  run on 30,000 rows); all three finished and all three were used.
- **Four offline experiments, not submissions:** FRIGID de novo (50 molecules, 3819 generations of
  which 3588 valid, all 50 in the denominator), the FPNet encoder check (50 held-out molecules), the
  MIST `msg_aug` probe (50 FRIGID molecules, 0.385 against the 0.40 gate) and the external-library reach
  check (3,096 new structures, 0 / 400 test molecules reached); see above.
- **One change not run, on a reach check:** PubChem promotion at `rel > 100` on E006. Decided before
  any submission, from public evidence: its public sibling notebook already scored 0.444 against 0.446,
  so its best case sits inside the noise.
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
| E006 engine-2 fusion weight 0.8 (e006-v3) | **before** it was submitted | one setting per submission; scored 0.437, discarded |
| skip `rel > 100` on E006 | **before** any submission | reach check: its public sibling scored 0.444 against 0.446 |
| the spectra-fork's settings (depth-7 rankers, standardised candidates, version 53's seeds) | by its author, before we ran it | not ours |
| the 0.456 fork's settings (PubChemLite in engine 2, N1A, popularity 0.15) | by its author, before we ran it | not ours; the popularity weight was chosen on the public leaderboard |
| the data fork's settings (proxy re-ordering, fold-safe FPNets, sixth PubChem slot, candidate-aware promotion) | by its authors, before we ran it | not ours; its lineage was scored on the public leaderboard |
| the neural fork's settings (frozen 0.465 baseline, spectrum model, RRF constant 3) | by its author, before we ran it | not ours; its weight 0.25 was scored on the public leaderboard |
| neural re-rank weight 0.50 (neural-fork) | **before** it was submitted | one setting per submission; scored 0.444, discarded |
| the merged fork's settings (re-rank weights 0.25 / 0.75, fragment check, `[P0, V0]` leader merge) | by its authors, before we ran it | not ours; each change was scored on the public leaderboard (0.472 each) |
| the dedup fork's settings (prior weights 0.10 / 0.05 / 0.02, library gate 0.87, leader pin, InChIKey14 dedup) | by its authors, before we ran it | not ours; scored on the public leaderboard upstream (0.479) |
| seed blend: equal weights, shipped schedule, three seeds | **before** the retrains (2026-11-03) | one change per submission; scored 0.463, discarded |
| mass prior: V22's constants (centre −1.5 ppm, free 3, full 8, weight 0.15; rank-1 rule at 5 ppm, library < 0.9) | by its authors; adopted **before** our run (2026-11-08) | one change per submission; band checked on train first; scored 0.479, discarded |
| the extlib fork's settings (scorer ensemble 0.5 / 0.5, mass prior, documentation prior, extlib cosine 0.8, K = 4) | by its authors, before we ran it | not ours; scored on the public leaderboard upstream (0.487) |
| isomer reranker for v2: `V8 = 'BOTH'`, `RR8_ALPHA = 0.20` | **before** its run (2026-11-13) | one change per submission; small weight chosen so it can only reorder same-formula isomers; 0.487, no change; v3 lowers it to 0.15 (fixed before its run), running |
| MIST `msg_aug` gate: mean Tanimoto ≥ 0.40 on the 50 FRIGID molecules | **before** the run (2026-11-01) | pre-registered; 0.385 fails it, so no FRIGID decode was run |
| external-library reach: cosine ≥ 0.5, 20 ppm, same ion mode, structures absent from train | **before** the count | a reach check; 0 / 400, so no submission |
| analogue propagation (fork477 v3: anchor ≥ 0.87, cosine ≥ 0.70, λ 0.10, rank 1 frozen) | **before** it was submitted | its reach (113 / 400) checked first; no offline panel can score it; scored 0.477, no change, discarded |
| FRIGID go bar: top-10 exact ≥ 5 / 50 with the true formula | **before** the run (2026-10-18) | pre-registered; MIST-MSG over CANOPUS also fixed before |
| FRIGID encoder diagnostics (MassBank vs Enveda Tanimoto, oracle fingerprint) | **after** the 0 / 50 result | post hoc; explain the failure, not a statistic |
| FPNet check gates: 4096-bit Morgan space, mean Tanimoto ≥ 0.45 | in the check's protocol, before its results | sample re-drawn after finding only 10 / 50 held out; threshold sweeps, oracle top-k, trained-vs-held-out and 50-way retrieval are post hoc |
| keep the attached GLACIER checkpoint (ms-pred's public MassSpecGym-trained) | 2026-10-17, before any swap | the paper's contrastive fine-tune gives 69.95 against 69.69 MassSpecGym top-1; its NIST gain needs NIST-licensed weights |

Choices marked **after** were selected on the public leaderboard, so our pipeline's 0.356 is an
optimistic estimate of the private score; the forks' 0.415, 0.446, 0.456, 0.465, 0.477 and 0.479 involve none of ours, though the E006 threshold and the 0.456 notebook's popularity weight were tuned on the public split by their authors. Differences under about ±0.005 — the run-to-run noise
of identical runs — were treated as no change.

## What is not measured

- **The private score.** Every number here is public-split.
- **Why fork-v1 scores 0.415 on our account against 0.421 public.** Both our runs agree, so the
  0.006 is not run-to-run noise; the E006 fork reproduced its public 0.446 exactly on the same account
  and day, so it is not a general account or day effect. The cause is not established.
- **Why the 0.451 notebook scores 0.445 for us.** Its 0.451 is for version 52; we could fork only
  version 53, which adds ranker seeds. Seeds or non-reproduction: not established.
- **Whether a de novo route could help with a better encoder.** FRIGID found 0 / 50 with MIST's
  fingerprints; an encoder fine-tuned on Enveda's spectra is untested. The 39–52% of the test outside
  PubChem and COCONUT is the SOTA sweep's estimate, not a measurement on the hidden test.
- **Whether the 0.479 layer helps at all.** +0.002 over 0.477 is inside the source's ~0.002 run-to-run spread; the
  layer's effect on its own (documentation prior, leader pin, dedup) is not measured by us.
- **Why the seed blend lost.** Untested: the seeds agree closely with the public scorer and two of three validate
  lower, so equal weights may just dilute it; other weights, or seeds from a different recipe, were not run.
- **What in the 0.487 notebook earns its +0.008.** We ran it unchanged; the split between the scorer ensemble, the
  external-library promotion and the other layers is not measured, and +0.008 itself is inside the noise. How many
  hidden-test molecules the 324,884 external spectra promote is not visible: the smoke run's molecules are train copies.
- **How many hidden-test answers the mass prior reached.** The commit run scores only the 12 smoke molecules; its
  effect on the public split is a measured zero, but the split between the ranks 2–25 penalty and the rank-1 rule is not.
- **What `msg_aug` would do to ranking.** At Tanimoto 0.26 FRIGID found 0 / 50; whether 0.385 gives any hits, or helps
  as an FPNet replacement or a re-ranking feature, was not run.
- **What in the 0.477 notebook earns its +0.012.** We ran it unchanged; the split between the leader-safe
  re-rank and the fragment check (0.472 each upstream) is unmeasured by us.
- **Whether a larger λ would help analogue propagation.** At λ 0.10 it gave no change (0.477 against 0.477; at most 113 / 400 molecules reachable); a larger λ is untested, and why it changed nothing (few lists moved, or the moves did not matter) is not established.
- **What in the version-26 notebook earns its +0.009.** We ran it unchanged; the split between the proxy
  re-ordering, the fold-safe FPNets, the sixth PubChem slot and the candidate-aware promotion is unmeasured by us.
- **What in the 0.456 notebook earns its +0.010.** We ran it unchanged; the split between PubChemLite,
  N1A and the popularity weight rests on its author's public scores.
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
  copied spectra were built deliberately, not by library duplication. A public labelled panel does not
  calibrate either: sanpier's 120-molecule panel ranks the 0.477 recipe below the 0.465 baseline
  (−0.004 ± 0.003), where the leaderboard gave +0.012, so the leaderboard remains the only measurement.
- **Why the data is shaped as it is.** 1189 / 1213 test spectra are bit-identical to a train row,
  covering every molecule, yet that row's structure is the answer for at most ~46 / 400.
- **Candidate recall.** Whether the answer is among the ~1038 candidates per molecule is unknown;
  rank and recall cannot be separated from one MRR.
- **29 / 400 molecules get no forward-model evidence in our pipeline** (no positive `[M+H]+` / `[M+Na]+`
  spectrum); the forks' coverage of the hidden test is unknown.
- **The image and RDKit switch** (v11 → v12) was not measured on its own; v11's output was lost to
  a session restart before it could be compared.
