# CASMI: library search + learned ranker — design

Status 2026-10-13: plan, nothing built. Follows the 0.421 comparison in
[PR #85](https://github.com/charleneleong-ai/ai4science/pull/85).

## Why

Our best kept pipeline scores 0.356; the public notebook our components come from scores 0.421. The
comparison of the two found the largest difference in the base ranking:

- **0.421** searches every test spectrum against all of train, propagates library similarity to
  structurally similar candidates, and ranks a pool of train ∪ COCONUT ∪ generated derivatives with a
  public LightGBM ranker (`ranker_v4b.pkl`, 4 boosters, 160 features). PubChem candidates are injected
  only when the library match is weak (`lib_max < 0.7`), and the forward models are switched off for
  library hits.
- **Ours** excludes the train structure behind each copied spectrum and ranks z(FPNet f·z) + 0.25 ×
  popularity over a popularity-capped PubChem pool.

Our exclusion rests on a visible-test statistic: the visible test is train examples (the 0.421
notebook's own comment, `src.py:133`), while the leaderboard scores a hidden test. That is why every
copy-targeted change measured as noise. On the hidden test, library matches are evidence to use.

## Options

**A. Rebase on the 0.421 notebook.** Fork it unchanged (all 15 datasets attached, T4), submit it once to
reproduce ~0.421 under our account, then make one change per submission on top of it.
- One submission to land about +0.065, if it reproduces.
- Our own code becomes the ablation record; later changes are measured against 0.421, not 0.356.
- Costs: a 161 KB notebook we did not write; its runtime on the full hidden test is unknown to us
  (the 9 h limit applies).

**B. Port the base ranker into our notebook.** Attach `casmi26-v4b-models` and `casmi26-v2-pool`, build
the `Engine` / `V1FE` / `rank_score` stack as `src.py:304-348` does, and replace our candidate pool and
FPNet + popularity score with its base lists (top 60, key-deduplicated); keep our second PubChem pass,
ICEBERG and GLACIER on top. Drop the copy exclusion and the similarity filter in the same step, since the
engine does not exclude.
- Keeps our pipeline, tests and gates; the notebook stops being self-contained (it imports the v4b
  package from an attached dataset).
- Then add, one submission each: PubChem injection gated on `lib_max < 0.7` (`pc_join.py`); forward
  models off for library hits (`fusion_core.py:144-148`); fragment rescoring for molecules ICEBERG
  cannot score (`frag_rescore.py`); engine-2 RRF fusion last, as the heaviest.
- Five or more submissions to approach 0.421, each attributable.

## Recommendation

**A first, as a measurement.** One submission tells us whether 0.421 reproduces on our account and how
long the full run takes. If it does, building on it is the fastest route to the leaderboard top (0.48),
and B's steps become ablations of 0.421 rather than reconstructions of it. If it does not reproduce,
B is the fallback, and the failed reproduction tells us something about the hidden test.

## Gate for each submission

As now: the model and ranker loaded (`engine ready; ranker features 160 boosters 4` in the log), the
base lists have no errors (`base lists done, errors 0`), ICEBERG and GLACIER returned scores, output is
400 × 25. G1 is unchanged: beat the best score so far by more than the measured ±0.003 noise.

## Not established

- Whether the 0.421 score reproduces on a different account and day (hidden-test reruns may vary).
- Runtime of the full scoring rerun; the commit run only scores 12 smoke molecules.
- How many hidden-test molecules are library hits; the visible test cannot tell us.
- Both options need the weekly GPU quota, which is exhausted until it resets.
