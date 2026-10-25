# fork466

A fork of the public notebook
[PubChem + Neural025 — LB 0.466](https://www.kaggle.com/code/sanpier/pubchem-neural025-lb-0-466) by Saner Turhaner,
run as the private kernel `charyeezy/casmi-neural-fork`. `neural_fork.ipynb` was vendored byte-identical to its
current version as of 2026-10-25 (md5 `af68cfeeb06a251186cf043684ca3be2`); [`kernel-metadata.json`](kernel-metadata.json)
sets our kernel id, private, and drops an empty dataset entry.

The notebook runs a frozen 0.465 baseline once (embedded as base64, sha256-checked), then re-orders each molecule's
final 25 by reciprocal-rank fusion of the baseline rank with the rank from a separately trained spectrum model
(`sanpier/casmi-spectrum-model-v1`): `1/(3 + baseline_rank) + w/(3 + neural_rank)`. Candidate sets never change.
Its frozen baseline is the 0.464 lineage with the candidate-aware rank-1 PubChem promotion but one fold-safe FPNet
and no sixth PubChem slot; our [`fork464`](../fork464) (v26, five folds and the sixth slot) scores the same 0.465.

| upstream version | neural weight `w` | public MRR@25 |
|---|---|---|
| reported | 0.25 | 0.466 |
| current (vendored) | 0.15 | not scored upstream |

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | `NEURAL_WEIGHT` 0.15 → 0.50 (cells 1 and 9; the report's experiment name follows) | 0.444 (submission 57615727): −0.022 against weight 0.25 upstream and −0.021 against our 0.465; discarded |

Run with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh); gate `~/casmi-data/gates/gate_neural.py` (the
0.451-lineage checks on the frozen baseline, plus the neural report at weight 0.50 and its audit showing the smoke
molecules processed and re-ordered).
