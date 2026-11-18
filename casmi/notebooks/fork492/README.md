# fork492

A fork of alex chilton's public notebook
[CASMI26 0.492 | What Moved the Board: 59 Arms](https://www.kaggle.com/code/alexchilton/casmi26-0-492-what-moved-the-board-59-arms)
(version 2), run as the private kernel `charyeezy/casmi-fork492`. The pipeline is [fork487](../fork487)'s, plus
one cell, N-1c: a test-set molecular-network rerank. Test spectra with modified cosine ≥ 0.6 on ≥ 6 peaks are
linked, and each list is reordered by z(−log rank) + `LAM` × z(neighbour score), with `LAM` raised from 0.25 to 0.35.
The source page verifies **0.492** (submission 58380172). Its own LAM sweep scored 0.487 / 0.487 / 0.488 / 0.492 /
0.484 at 0.0 / 0.15 / 0.25 / 0.35 / 0.50. So the 0.492 is the best of five, about one molecule above fork487.

The one edit to [`fork492.ipynb`](fork492.ipynb) is `COMMIT_MODE = 'fast'` → `'smoke'`, 6 bytes. In `'fast'` the
commit run skips the pipeline, leaving nothing to gate. [`kernel-metadata.json`](kernel-metadata.json) has our
private kernel id, with `id_no` and `keywords` dropped.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | `COMMIT_MODE` smoke only | **0.492** (submission 58506761, 2026-11-18; reproduces the source; final pick A) |

The gate is `~/casmi-data/gates/gate_fork492.py`: `gate_fork487.py`, then the `[n1c] LAM=0.35` line and no
`[n1c] FAILED`, because N-1c leaves the submission unchanged on any exception.

Launched as the fallback after [fork495h](../fork495h) could not be submitted. The commit run's N-1c line was
`[n1c] LAM=0.35 TAU=0.6 molecules 400, with spectra 391, with neighbors 114, edges 109, changed lists 2, top1 changed
0, 80s`. The +0.005 over fork487's 0.487 is about one molecule, inside the ±0.016 noise. It is promoted as the new
base, not counted as a measured gain.
