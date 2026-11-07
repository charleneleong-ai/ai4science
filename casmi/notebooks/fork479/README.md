# fork479

An unchanged fork of the public notebook
[CASMI26 LB 0.479 | MIDX Exact Dedup Fill](https://www.kaggle.com/code/alexchilton/casmi26-lb-0-479-midx-exact-dedup-fill)
by alex chilton, which itself runs Baidalin Adilzhan's *casmi26 midx exact dedup fill* unchanged (that page was no
longer reachable through the API on 2026-11-01). It runs as the private kernel `charyeezy/casmi-fork479`.
[`fork479.ipynb`](fork479.ipynb) is byte-identical to the source as of 2026-11-01 (md5
`a6871ca326db69f8c4a5caff8d6b7d36`); only [`kernel-metadata.json`](kernel-metadata.json) differs (our kernel id,
private). It attaches the source's eight public datasets: the seven of [`fork477`](../fork477) plus
`metric/rdkit-2026-3-3-wheel`. It also uses the source's GPU image.

It is the [`fork477`](../fork477) pipeline (0.477) with a final layer on ranks 2–25:

- **Documentation prior:** PubChem popularity 0.10, a spectral-library bonus 0.05 (only where the best library
  match reaches 0.87), and a natural-product bonus 0.02, all read from `pool_popularity.csv`.
- **The old leader pinned at rank 2** after the fragment check moves rank 1.
- **Metric-exact InChIKey14 deduplication:** stereoisomers that the scorer counts as the same answer give their slots
  to new skeletons, and empty slots are filled with scored outside candidates.

The source reports 0.479 on its author's account. Its own table shows the same code scoring about 0.002 apart between
runs, and two notebooks reported as 0.478 on their pages scored 0.468 and 0.476 when rerun. So +0.002 over 0.477 is
inside the noise.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: upstream as of 2026-11-01 | **0.479** (ref 57867559, 2026-11-01; reproduces upstream; KEEP, +0.002 over 0.477 is inside the noise) |
| 2 | ranks 2–25 neural score = mean per-molecule z of the public scorer and three seed retrains of it ([`b9627b0`](https://github.com/charleneleong-ai/ai4science/commit/b9627b0)); private dataset `charyeezy/casmi26-seed-scorers` | 0.463 (ref 58046573, 2026-11-06; DISCARD, −0.016) |

Run with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh). The gate is `~/casmi-data/gates/gate_fork479.py`:
the `fork477` gate, plus the final layer's own log lines, and no `Note: ... exception` line. Those lines are how the
layer's fallbacks show: if `pool_popularity.csv` or the neural evidence fails to load, it proceeds with a neutral
prior and does not fail the run.
