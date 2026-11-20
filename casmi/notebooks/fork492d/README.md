# fork492d

[fork492](../fork492) (public 0.492: fork487 plus the test-set network rerank at `LAM` 0.35) plus one cell: DreaMS
on ranks 2–25. Run as the private kernel `charyeezy/casmi-fork492d`.

The new cell (code cell 16 of [`fork492d.ipynb`](fork492d.ipynb), headed `DreaMS on ranks 2-25`) runs after the APEX
library overlay, so after N-1c.
- **Spectrum side:** the official frozen DreaMS embedding model embeds each positive-mode test spectrum. It uses
  DreaMS's own preprocessing: the top 100 peaks in m/z order, base-peak relative, with the precursor first at
  intensity 1.1. The cell averages the unit embeddings over the molecule's spectra. Molecules with only negative
  spectra are skipped.
- **Molecule side:** a tower maps each candidate's 4,096-bit Morgan r2 fingerprint into the same 1,024-d space.
- **Reorder:** ranks 2–25 by 0.5·z(−log rank) + 0.5·z(cosine). Rank 1 never moves, and the cell asserts it. On any
  error the submission is unchanged.
- **Tower:** [`dreams_tower.py`](../../scripts/dreams/dreams_tower.py), run on the A100 by
  [`run_dreams.sh`](../../scripts/dreams/run_dreams.sh).
  - An MLP (4,096 → 2,048 → 2,048 → 1,024) trained with symmetric InfoNCE against each molecule's mean DreaMS
    embedding (positive mode, at most 4 spectra, timsTOF first).
  - Trained on `train.parquet` only: the natural-product holdout's `train` split, 250,339 molecules. Early stopping
    on its `validation` split, 957 molecules: best retrieval MRR 0.571 over 30 epochs. The `eval` split and its near
    neighbours are excluded.
- **Dataset:** the private, versioned dataset `charyeezy/casmi-dreams-frozen` holds
  `DreaMS_embedding_model_torchscript.pt` and its settings from
  [HF roman-bushuiev/DreaMS](https://huggingface.co/roman-bushuiev/DreaMS) (MIT, unchanged, used frozen), plus
  `tower.pt`.

**Licence and ruling.** DreaMS is MIT. The host ruled in
[discussion 754610](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754610):
- 5 Nov: "Officially released DreaMS checkpoints are fine to use."
- 17 Nov: `ssl_model.ckpt` also counts as an official checkpoint. That answered a question describing exactly this
  use: frozen encoder weights shared as a Kaggle dataset, with a small head trained only on the competition data.

The DreaMS weights are not fine-tuned. Only the tower is trained.

**Reach,** replayed on the 953 natural-product holdout lists, with the shipped scorer's order standing in for the
pipeline:
- 915 lists have a DreaMS embedding, and 839 are reordered;
- the true structure is in ranks 2–25 in 477 lists and moved in 259, 198 up and 61 down;
- MRR@25 0.576 → 0.595, +0.019 (SE 0.002): Enveda-180 +0.013 (727 lists), natural products +0.038 (226).

That holdout does not predict the leaderboard, so this shows reach, not a gain.

The gate is `~/casmi-data/gates/gate_fork492d.py`: `gate_fork492.py`, then DreaMS's `[dreams] done:` line with rank 1
unchanged, nulls at most 5% of molecules, and at least one list changed, and no `[dreams] FAILED`.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | fork492 + DreaMS tower cosine on ranks 2–25 (weight 0.5) | **0.486** (58580239), −0.006: DISCARD |

**Result (2026-11-19).** The commit run's gate passed:
`[dreams] done: {"molecules": 400, "spectra": 987, "no_positive": 27, "short": 362, "scored": 11, "null": 0, "changed": 11, "top1_changed": 0, "seconds": 19.3}`.
Submission 58580239 scored 0.486 against fork492's 0.492, inside the ±0.016 noise but lower, so it is discarded. It
is the same score as MIST's ([`fork492m`](../fork492m)), whose offline gain was twice as large. Lower weights,
negative-mode spectra and fine-tuning DreaMS were not tried.
