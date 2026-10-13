# fork421

An unchanged fork of the public notebook
[CASMI26 | v4n Fusion + PubChem on public 0.421](https://www.kaggle.com/code/huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421)
(public score 0.421), run as the private kernel `charyeezy/casmi-v4n-fork`. `v4n_fork.ipynb` is
byte-identical to that notebook as of 2026-10-13; only [`kernel-metadata.json`](kernel-metadata.json) differs
(our kernel id, private). Its purpose is option A of the
[library search + ranker design](../../docs/specs/2026-10-13-casmi-library-ranker-design.md): reproduce 0.421 on
our account, then change one thing per submission on top of it.

A commit run scores only 12 smoke molecules (the visible test is train examples) and writes `CCO` for the rest;
the scoring rerun processes the hidden test in full.
