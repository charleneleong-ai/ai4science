# fork446

An unchanged fork of the public notebook
[CASMI26 | E006 + PubChem Promo rel>150 | LB 0.446](https://www.kaggle.com/code/seyitkaangunes/casmi26-e006-pubchem-promo-rel-150-lb-0-446)
(public score 0.446), run as the private kernel `charyeezy/casmi-e006-fork`. `e006_fork.ipynb` is
byte-identical to that notebook as of 2026-10-17 (md5 `9829a9667030e7044027f0cd33fd7f5d`); only
[`kernel-metadata.json`](kernel-metadata.json) differs: our kernel id, private, and the GPU docker image `fork421` runs on. The source's
metadata names Kaggle's CPU image (`gcr.io/kaggle-images/python`), on which the notebook's T4 preflight fails. It descends from the same
lineage as [`fork421`](../fork421): flexonafft's v80 recipe (0.436) → E006, which fuses a second engine's
ranks and adds ChEBI/LipidMaps candidates and simulated ranker rows (0.442) → PubChem promotion at
`rel > 150` instead of 250 (0.446). It attaches 16 public datasets, 14 of them shared with `fork421`.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none, on the source's CPU image | failed at preflight (`nvidia-smi` exit 12, no GPU); not submitted |
| 2 | none: reproduction of 0.446 on the GPU image | **0.446** (submission 57298936): reproduces exactly |
| 3 | engine-2 weight in the rank fusion `ALPHA` 0.6 → 0.8 (cell 35; every fusion call reads it) | 0.437 (submission 57310111): worse, discarded |

The run gate is `~/casmi-data/gates/gate_fork446.py`: the engine and ranker loaded, the engine-2 cell
applied, base inference had 0 errors, no stage logged a failure, and the output is 400 rows with the 12
smoke molecules full. A commit run scores only those 12 (the visible test is train examples); the scoring
rerun processes the hidden test in full.
