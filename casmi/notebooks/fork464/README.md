# fork464

An unchanged fork of the public notebook
[CASMI 26: The Data Behind Spectra to Structures](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-the-data-behind-spectra-to-structures)
by dmitriigluzdov, run as the private kernel `charyeezy/casmi-data-fork`. `data_fork.ipynb` is byte-identical to its
current version 26 as of 2026-10-24 (md5 `4edf7f5dbb3bee2aff0cc52f63961388`); only
[`kernel-metadata.json`](kernel-metadata.json) differs (our kernel id, private). It attaches five of the author's
datasets (the fold-safe FPNet among them) and the GPU image the other forks use. The notebook's own version table,
from the method of our [`fork456`](../fork456):

| upstream version | change | public MRR@25 |
|---|---|---|
| 13 | GLACIER proxies for molecules without an [M+H]+ spectrum | 0.456 |
| 16 | proxies kept where ICEBERG scores exist; joint ICEBERG and GLACIER re-ordering of proxy-covered lists | 0.464 |
| 19–20 | five fold-safe FPNets in the PubChem channel; a sixth PubChem slot at rank 24 | 0.464 |
| 26 | rank-1 PubChem promotion judged on the leading proposal's own scores ([Saner Turhaner](https://www.kaggle.com/code/sanpier/candidate-aware-pubchem-promotion-0-465), 0.465 there) | not scored upstream |

Kaggle serves only the current version, so this reproduces version 26, not a scored version.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: upstream version 26 | running |

Run with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh); gate `~/casmi-data/gates/gate_fork451.py`.
