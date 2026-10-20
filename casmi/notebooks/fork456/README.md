# fork456

An unchanged fork of the public notebook
[casmi26-v465-v15-pop015](https://www.kaggle.com/code/guhongbin/casmi26-v465-v15-pop015) by Hongbin Gu (public score
0.456), run as the private kernel `charyeezy/casmi-pubchemlite-fork`. `pubchemlite_fork.ipynb` is byte-identical to it
as of 2026-10-20 (md5 `9d1ab164357187405aa67fd1335d081d`); its code cells also match line for line the documented copy
[CASMI26 | LB 0.456 PubChemLite + N1A, explained](https://www.kaggle.com/code/omerfarukmerey/casmi26-lb-0-456-pubchemlite-n1a-explained),
whose table gives the steps from the 0.451 notebooks:

| step | change | public MRR@25 |
|---|---|---|
| base | E006 + rel 150 + rank-2 hedge + standardised ChEBI / LIPID MAPS (cf. [`fork451`](../fork451)) | 0.451 |
| 1 | PubChemLite structures in engine 2's library | 0.455 |
| 2 | N1A: GLACIER scores for molecules without an [M+H]+ spectrum, from their other adducts | — |
| 3 | PubChem popularity weight 0.25 → 0.15 | 0.456 |

Only [`kernel-metadata.json`](kernel-metadata.json) differs (our kernel id, private); it attaches the same 15 public
datasets and GPU image.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: reproduction of 0.456 | running |

The run gate is `~/casmi-data/gates/gate_fork451.py` (the ranker loaded, the merge self-test passed, base inference had
0 errors, no stage logged a failure, 400 rows with the 12 smoke molecules full).
