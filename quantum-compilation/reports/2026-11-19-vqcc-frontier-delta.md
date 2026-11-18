# VQCC frontier delta, 2026-11-12 to 2026-11-19

Delta against [`2026-11-12-vqcc-frontier-delta.md`](2026-11-12-vqcc-frontier-delta.md). The
leaderboard API was unauthenticated, so board figures are manwithacat's briefing of 2026-11-18
([760157](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/760157)),
third-party.

## What changed

**The public frontier went from 25.405 to 28.21 (+11%) in a week, mostly in invsqrt.** Top 31.617
(was 29.52), gold line 28.753 (26.80), median 25.150, so our 25.366 had dropped below the median.
We have since submitted the public best as v4: leaderboard **28.18862**
([`58528892-parts.yaml`](../docs/experiments/58528892-parts.yaml)).

| task | v3 (ours) | best public, Apache-2.0 | gain |
| --- | ---: | --- | ---: |
| invsqrt | 22.13 (Q 62) | **27.91**, `d0227f3c`, Q 44, cost 74,405.67 ([maoxin123](https://www.kaggle.com/code/maoxin123/qcc-28-20-4-022-toffoli-chemistry-loaders)) | ×1.261 |
| dqi | 65.41 (Q 892) | **70.79**, `e208d350`, Q 890, 3.4525e11 ([Ertuğrul Özer](https://www.kaggle.com/code/ozertuu/quantum-ai-fast-runner-starter-notebook)) | ×1.082 |
| dfthcsos | 11.28 | **11.39**, `6a094aa9`, (8, 143, 44), T 14,706, λ_eff 46.92 | ×1.010 |

Ratios are our local verifier's on the files themselves. Unpublished: Syed Asad Ali reports
dfthcsos 293,957,864 (about 11.42, loaders 4,007 per pass;
[760029](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/760029),
author-reported).

## Disclosures

- **Chemistry is near exhausted** (760029). GrantPrater bounds λ_eff ≥ √(559.13 × 3.914898) = 46.786
  for any shape or factorisation, matching our `r1` floor; the public sits 0.26% above it.
  Selvaraj Sengottiyan reports five negative levers (grouping, Givens tree, 9-bit angles, J = 4–6,
  Λ refits). Syed's best recipe gains 0.28% and he concludes more needs a new loader or factorisation.
  Apih Pih: a k10 loader re-plan reaches 190 Toffolis against 189 today
  ([759820](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759820)).
- **DQI moves by merging patch sets** ([759760](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759760),
  author-reported). Open items each worth 0.1–0.2% of the task: 686 closing Toffolis if a second idle
  helper wire exists, a guard for the zero-word input class behind ~0.33 lost shots per draw, and
  ~7 Toffolis per shot from commuting control products into measurement-conditioned CZ/CCZ.
  The key-equation stage is reported at ~87.6% of Toffolis.
- The loader eligibility question (757522) and milestone timing (760066) have no host answer.
  Milestone snapshot 2026-11-20 23:59 UTC; final 2027-01-20 23:59 UTC.

## Papers

Nothing new on THC loaders or DQI decoders. [arXiv:2605.20334](https://arxiv.org/abs/2605.20334)
(*Halving the cost of QROM*) folds each unload into the next load; unread, and public loaders
already load XOR deltas between batches, so its novelty for us is unchecked.

## What it changes in the plan

- **invsqrt is not closed as a task.** Only our lookup-redesign route is (`r7`); the public Q 44
  line is ×1.26 past our old target by another route, and is now the baseline.
- **Chemistry stays closed** to tuning; only a new loader architecture could move it.
- **DQI has the most concrete open items**, each small.
- Not established: how the leader's ×1.12 over public splits across the three tasks.
