# fork421

An unchanged fork of the public notebook
[CASMI26 | v4n Fusion + PubChem on public 0.421](https://www.kaggle.com/code/huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421)
(public score 0.421), run as the private kernel `charyeezy/casmi-v4n-fork`. `v4n_fork.ipynb` is
byte-identical to that notebook as of 2026-10-13 in its first commit (kernel version 1); only
[`kernel-metadata.json`](kernel-metadata.json) differs (our kernel id, private). Each later version changes one
setting:

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: reproduction of 0.421 | 0.415 (submission 57288940) |
| 2 | `PC_N1=10000` (PubChem channel re-scores the top 10,000 pass-1 rows, not 5000; probe v3: 255 vs 236 / 280 held-out answers) | scoring |

Its purpose is option A of the
[library search + ranker design](../../docs/specs/2026-10-13-casmi-library-ranker-design.md): reproduce 0.421 on
our account, then change one thing per submission on top of it.

A commit run scores only 12 smoke molecules (the visible test is train examples) and writes `CCO` for the rest;
the scoring rerun processes the hidden test in full.
