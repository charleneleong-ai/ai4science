# fork477

An unchanged fork of the public notebook
[LB 0.477 · Two 0.472 Fixes, Merged the Right Way](https://www.kaggle.com/code/busyaprime/lb-0-477-two-0-472-fixes-merged-the-right-way)
by Akmal Xodarev (Busya PRIME), run as the private kernel `charyeezy/casmi-merged-fork`. `merged_fork.ipynb` is
(kernel version 1) byte-identical to its current version as of 2026-10-28 (md5 `74efef20daa877c43243024d830ef032`), the version its
page reports as submission 57722897 at public 0.477; only [`kernel-metadata.json`](kernel-metadata.json) differs (our
kernel id, private; the source's own GPU image, newer than the one our other forks pin). It attaches the same seven public datasets as [`fork466`](../fork466).

It runs the pinned 0.465 pipeline once, then stacks two public changes that each showed 0.472:

- **Leader-safe neural re-rank** (alex chilton): the spectrum model of [`fork466`](../fork466) at weight 0.25 for
  rank 1 and 0.75 for ranks 2–25. Our weight-0.50 run that let it move rank 1 everywhere scored 0.444.
- **Fragment check** (Roman Tamrazov's Spectral Parallax): an untrained MetFrag-style re-scorer that may replace the
  leader. Its leader moves are taken from the 0.50 (V2) order it was built on and laid over the 0.75 order.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | none: upstream as of 2026-10-28 | scoring (ref 57822991) |
| 2 (our v3) | one new cell after the merge: analogue propagation across test molecules | running |

**v3, analogue propagation.** The host allowed automated similarity between test molecules on 2026-10-26
(discussion 752297). Test molecules whose leader has a library match ≥ 0.87 are anchors. For every other molecule
whose spectra reach cosine ≥ 0.70 with an anchor, candidates below rank 1 gain 0.10 × cosine × Tanimoto(candidate,
anchor leader), and the tail is re-sorted by 1/(3 + position) + bonus. Rank 1 never moves (ranks 1–2 where the merge
moved the leader), and the candidate set never changes. Reach on the visible test: 113 of 400 molecules have a
neighbour at cosine ≥ 0.70, so that many at most can change. On `enveda-np-examples` such neighbours share structure
(Tanimoto ≥ 0.5) 46% of the time, against 8% below cosine 0.5; on `enveda-180` they do not (3%). No offline panel
can score it: sanpier's 120-molecule labelled panel ranks the 0.477 recipe below the 0.465 baseline, so the leaderboard
is the only measurement, and the expected change is inside its ±0.010 noise. In the commit run the 12 smoke molecules
have no neighbour, so the stage applies but changes nothing; the gate checks it ran and that no leader moved.

Run with [`scripts/kaggle_chain.sh`](../../scripts/kaggle_chain.sh); gate `~/casmi-data/gates/gate_fork477.py` (the
0.465 pipeline's checks read from its plain-text stage logs, both layer reports written, the fragment check not
skipped, and 400 rows with the 12 smoke molecules full); v3 adds `gate_fork477_v3.py` (the stage applied, and no
leader or candidate set changed against `submission_pre_analogue.csv`).
