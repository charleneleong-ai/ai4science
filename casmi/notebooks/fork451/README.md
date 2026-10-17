# fork451

An unchanged fork of the public notebook
[CASMI 26: From Spectra to Structures](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-from-spectra-to-structures)
by dmitriigluzdov, run as the private kernel `charyeezy/casmi-spectra-fork`. `spectra_fork.ipynb` is
byte-identical to that notebook's current version 53 as of 2026-10-17 (md5 `951a34f88affa3b89d7fdefc166e60f5`); only
[`kernel-metadata.json`](kernel-metadata.json) differs (our kernel id, private). It is the E006 pipeline of
[`fork446`](../fork446) plus three later changes, per the notebook's own version table (section 14):

| upstream version | change | public MRR@25 |
|---|---|---|
| 49 | PubChem promotion at rel > 150 (= `fork446`) | 0.446 |
| 51 | depth-7 trees in the second-engine rankers | 0.450 |
| 52 | standardised ChEBI and LIPID MAPS candidates in the second engine | 0.451 |
| 53 | eight and six random seeds for the second-engine rankers | not scored upstream |

Kaggle serves only the current version, so this fork reproduces version 53, not the scored version 52. Its
four datasets are the author's own repackaging of the inputs `fork446` attaches.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: upstream version 53 | running |

The run gate is `~/casmi-data/gates/gate_fork451.py`: the ranker loaded, the merge self-test passed, base
inference had 0 errors, no stage logged a failure, and the output is 400 rows with the 12 smoke molecules full.
