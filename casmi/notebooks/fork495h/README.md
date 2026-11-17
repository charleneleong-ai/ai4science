# fork495h

An unchanged fork of Heliosli's public notebook
[casmi26-0.495-public demo](https://www.kaggle.com/code/heliosli/casmi26-0-495-public-demo) (version 2, 53 votes on
2026-11-17), run as the private kernel `charyeezy/casmi-fork495h`. [`fork495h.ipynb`](fork495h.ipynb) is
byte-identical to the source (md5 `c17a664085795ae81a6ee72b7b444939`). Only [`kernel-metadata.json`](kernel-metadata.json)
differs: our kernel id, private, with the source's `id_no` and `keywords` dropped (a copied `id_no` makes the push
fail with a 403). It already ships `COMMIT_MODE = 'smoke'`.

**Why fork it:** it is the only public pipeline outside the fork487 lineage whose inputs are all public. Heliosli's
seven datasets went public on 2026-11-16, and no kernel source is needed. That makes it a candidate diverse second
pick. **The 0.495 is in the title only.** The page text is busyaprime's 0.477 write-up and cites no 0.495
submission, so the claim is unverified.

What it runs: the 0.477 merged pipeline ([`fork477`](../fork477)), then Heliosli's overlays, each of which leaves the
submission unchanged on failure:
- FRIGID + E1 de novo candidates on the second T4, merged after the rankers;
- RR3, a same-formula isomer reranker with a training-structure correction;
- a p1 gate and CV-JOINT, a joint ranker that re-judges rank 1;
- M5 parts A and B: an in-list isomer judgement, here in `swap` mode. Part B's `m5d1` methane-probe mode is off.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: upstream version 2 as of 2026-11-17 | running |

The gate is `~/casmi-data/gates/gate_fork495h.py`:
- `gate_fork477.py` for the base, with 400 rows and the smoke molecules full;
- each overlay's own `done` line;
- no fallback line (`overlay failed`, `self-test failed`, `skipped: …`, `all_methane`);
- M5 part B in `swap` mode.

Run with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh).
