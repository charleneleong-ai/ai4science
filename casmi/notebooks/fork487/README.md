# fork487

A fork of the public notebook
[casmi26-v589-alex0484-extlib](https://www.kaggle.com/code/guhongbin/casmi26-v589-alex0484-extlib) by Hongbin Gu,
whose code huichieh ran unchanged for a public 0.487 (submission 58185447, 2026-11-10). It runs as the private kernel
`charyeezy/casmi-fork487`. [`fork487.ipynb`](fork487.ipynb) is the source as of 2026-11-12 (md5
`ec45cc46672d00245bdb6bb0c39712f7`) with one change: `COMMIT_MODE = 'fast'` became `'smoke'`. In smoke mode the
commit run executes the whole pipeline on the 12 visible molecules, so the gate can see every stage. The scoring rerun
runs in full either way. [`kernel-metadata.json`](kernel-metadata.json) keeps the source's ten datasets, popwenman's
public `casmi26-scorer-train-v80` kernel and its GPU image, with our kernel id, private.

The pipeline is popwenman's Gold Apex V22 with four additions:

- **Two neural scorers:** the public spectrum model and an 80-epoch retrain ("OURS"), combined as
  0.5·z(public) + 0.5·z(ours) per molecule.
- **A timsTOF precursor-mass prior:** romanrozen's, as in [`fork479`](../fork479) v3.
- **A documentation prior:** PubChem, library and natural-product evidence on ranks 2–25, plus a LOTUS backfill.
- **External spectral libraries for rank-1 promotion:** Takumu Hata's extlib (200,724 spectra) plus the harmonised
  MassBank library, 324,884 reference spectra in all. A library match at cosine ≥ 0.8 can promote its structure to
  rank 1.

The +0.008 over 0.479 is one to two molecules on the roughly 124-molecule public split. That is inside the ±0.016
that two submissions differ by chance, so it is the new base (PROMOTE), not a measured gain over the 0.479 line.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none, apart from `COMMIT_MODE = 'smoke'` ([`2dc7d19`](https://github.com/charleneleong-ai/ai4science/commit/2dc7d19)) | **0.487** (ref 58272211, submitted 2026-11-12 23:06 AEDT; reproduces upstream; PROMOTE, +0.008 over 0.479, inside the noise) |
| 2 | v1 with the learned isomer reranker on: `V8 = 'BOTH'`, `RR8_ALPHA = 0.20` ([`28c5ce3`](https://github.com/charleneleong-ai/ai4science/commit/28c5ce3)); gate `gate_fork487_v2.py` | 0.487 (ref 58292050, submitted 2026-11-13 13:14 AEDT; DISCARD, no change; the commit run's reranker ran on 400 lists, skipped the 12 smoke lists as library matches and reranked 0, so its reach on the hidden test was not visible) |
| 3 | v2 with `RR8_ALPHA = 0.15` ([`62b95fb`](https://github.com/charleneleong-ai/ai4science/commit/62b95fb)); same gate | 0.487 (ref 58304138, submitted 2026-11-13 20:34 AEDT; DISCARD, no change; the commit run's reranker again ran on 400 lists and reranked 0. Three runs of this notebook, reranker off, α 0.20 and α 0.15, all scored 0.487, so the reranker line is closed) |

Run it with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh). The gate is
`~/casmi-data/gates/gate_fork487.py`, kept outside the repo.

- **What it requires:** the [`fork477`](../fork477) gate, plus these log lines:
  - `OURS BUNDLE READY`
  - the `ENSEMBLE neural = 0.5*z(public) + 0.5*z(ours)` line
  - `mass prior {"status": "on"`
  - `Documentation Prior: reordered`
  - the three `[apex-extlib]` lines `Loaded`, `Total references after Source B` and `SUCCESS`
- **What it rejects:**
  - `private model not found`
  - `Fast commit`
  - `overlay failed`
  - `loading error`
  - `off after an error`
  - any `Note: ... exception` line other than the LOTUS skip

On the v1 commit run, the extlib stage loaded 324,884 spectra and found 0 matches, as expected, since the smoke
molecules are train copies. The mass prior touched 394 of 400 lists, and the documentation prior reordered 12.

**Pushing the copy failed with a 403 at first.** `kaggle kernels push` answered "You don't have access… accept its
rules". The cause was the source's `kernel-metadata.json`, which carried its owner's `id_no`; that number overrides
`id`. Removing `id_no` (and `keywords`) fixed the push.
