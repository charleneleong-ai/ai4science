# VQCC frontier delta: where the leaders' advantage comes from (2026-11-12)

*Swept 2026-11-12, about 10:00–11:00 UTC. The raw excerpts and URLs are in
[`research_notes/2026-11-12-vqcc-frontier-delta/notes.md`](../research_notes/2026-11-12-vqcc-frontier-delta/notes.md).
Labels: **observed** (Kaggle CLI or Kaggle pages), **author-reported** (the poster's own claim), **ours** (our
arithmetic or inference).*

## Summary

- **The leaders have disclosed nothing technical.**
  - No prize has been awarded yet. The milestone snapshot is 2026-11-20 23:59 UTC; the single-task prizes are
    judged at the end.
  - No per-task leaderboard exists.
  - The rank-1 member's only forum post is about token cost.
- **The public frontier is 25.40482**, observed on the LB as submission 58260294. That is ×1.0015 on our
  25.36626.
- **The leader is now 29.52446**, Schrödinger's cat lived at 07:13 UTC. That is ×1.164 on our score, ×1.577 on the
  product of the three cost ratios.
- **The gap has been flat for 10 days.**
  - On 2 Nov the public frontier was 23.94 against a top score of 28.00, a gap of ×1.169.
  - Today it is 25.40 against 29.52, a gap of ×1.162.
  - Over those 10 days the public rose ×1.061 and the top ×1.055, so the leaders keep about one public fortnight
    ahead.
- **Two levers are closed.**
  - **Chemistry λ_eff**: we are within 0.2% of a proven floor.
  - **DQI width**: the 840 syndrome qubits are the *input* register, so 892 cannot fall below 840.
- **What remains is Toffoli and qubit compilation work.**
  - In invsqrt and/or the chemistry loaders: the leader's ×1.57 must come mostly from these (ours).
  - Possibly also DQI's Toffoli count.
- **There is one concrete public signal of private headroom.** A rank-31 player reports **invsqrt 25.1×**, against
  22.20 in public (author-reported, method undisclosed).

## Best known per-task ratios

| Task | Ours (58252192) | Best public (source) | Disclosed by anyone above public | Single-task ratio that alone closes the gap to 29.52 (ours) |
| --- | --- | --- | --- | --- |
| invsqrt | 22.1255 (93,863.4, Q 62) | **22.1984** (93,555.26, Q 62). File `9c708f7b`: asteriospanoras's splice of WCao24's normaliser and Thomas Gamet's prefix ([asteriospanoras](https://www.kaggle.com/code/asteriospanoras/lb-25-367-verified-portfolio-dqi-readout-fix), [attabeezy](https://www.kaggle.com/code/attabeezy/gq-portfolio-frontier)) | **25.1×**, mikelou1 ([756998](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756998), author-reported, no method) | 34.84 |
| dqi | 65.4149 (3.7361×10¹¹, Q 892) | **65.4275** (3.73535×10¹¹, Q 892, R 1.19126, about 394.1k Toffolis per shot). File `eae396ba`, Atakan Aldemir V15 ([iso-gap h5](https://www.kaggle.com/code/atakanaldemir/quantum-circuits-iso-gap-dqi-h5)) | not found | 102.70 |
| dfthcsos | 11.2771 (2.97558×10⁸) | **11.2893** (2.97237×10⁸; λ about 46.9, T 14,840, Q 427). maoxin123's 4,089-Toffoli loaders ([v22](https://www.kaggle.com/code/maoxin123/qcc-v22-4089-toffoli-chemistry-loaders)) with prvsiyan's V40 weights ([verifier-guided](https://www.kaggle.com/code/prvsiyan/verifier-guided-quantum-circuit-optimization)) | not found | 17.72 |
| **Blend** | **25.36626** | **25.40482** (LB, observed) | | **29.52446** (rank 1, observed) |

The best-public costs are **independently re-verified** with the unchanged official verifiers by Kirill Tymchishin
([retained-product audit](https://www.kaggle.com/code/kirilltym/quantum-circuits-retained-product-audit)). We have
not re-run them ourselves.

## What the top teams have disclosed

| Who | LB | What they said | Source |
| --- | --- | --- | --- |
| pjmathematician (competition_ranking 1) | 29.5 band | Uses Fable 5.2 for planning and Kimi for coding. No technique. | [757915](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757915) |
| GrantPrater (unranked, replying in the leader's thread) | — | "constraining the search space … [not] representing the problem in slightly the wrong way that leads to a combinatorial explosion" | [757915](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757915) |
| Some Joint (rank 5, 28.37) | 28.37 | The only top-10 team with a public notebook. Its [24.65 chemistry refit](https://www.kaggle.com/code/somejoint/quantum-circuits-24-65-local-gm-chem-refit) dates from Nov 4 and is far below their LB score. | kaggle CLI |
| mikelou1 (rank 31 on Nov 8) | about 26.5 | invsqrt **25.1×**, with public parts for the other two tasks | [756998](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756998) |

Check on the mikelou1 claim (ours): 25.1 with public DQI and chemistry blends to 26.46, which is consistent with rank
about 31.

Nothing in the forum, the prize pages or the notebooks gives a top team's per-task ratio. The question of **where**
the advantage sits can therefore only be answered by elimination.

**Elimination (ours).** The leader needs the product of the three cost ratios to be ×1.570 over the public
frontier.

- **Chemistry λ_eff: ≤ ×1.0021 on chemistry cost.** The bound λ_eff ≥ 46.79 holds for any factorization at the gap
  floor ([753633](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753633)).
  RezOps's ALS and LM refits found the fit floor to be a property of the table shape
  ([758129](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758129)).
- **DQI width: ≤ ×1.128 even at zero workspace.** The shipped verifier requires `S0..S69` (70 × 12 = 840 qubits) as
  the syndrome **input**, and every other qubit is workspace.
- **DQI's public line is nearly saturated.** It moved only ×1.011 between 2 Nov and 12 Nov (64.72 → 65.43), while
  invsqrt moved ×1.075 (20.65 → 22.20) and chemistry ×1.099 (10.27 → 11.29). Sources:
  [755038](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755038) and
  today's notebooks.
- **Conclusion.** Most of the ×1.57 must be Toffoli and qubit compilation in invsqrt and the chemistry block
  encoding, plus DQI's Toffoli count.
  - For example, if the leader's DQI is at the public level, invsqrt × chemistry must be ×1.57 over public. That is
    about +25% on each, e.g. invsqrt about 27.8 and chemistry about 14.1.
  - This is an inference, not a disclosure.

## Levers ranked by plausible gain on the blend

The blend factor is the cube root of the cost ratio on one task. A ×1.10 cost ratio gives ×1.032 on the blend,
×1.20 gives ×1.063, and ×1.25 gives ×1.077.

| # | Lever | Evidence that the headroom exists | Cost ratio (plausible) | Blend factor | Blend from 25.366 |
| --- | --- | --- | --- | --- | --- |
| 1 | **invsqrt: fewer Toffolis and qubits** (shifter CSWAPs, MBU on the phase-stage uncompute, fewer rotations) | mikelou1's 25.1× (author-reported). The public moved +7.5% in 10 days. The shifter is in the high 40s against a counting bound in the high 30s ([757030](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757030)). One qubit is 1.6% of cost ([753527](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753527)). Uncompute is about a quarter of Toffolis ([755302](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302), author-reported). | 25.1 / 22.1255 = ×1.134 | ×1.043 | **26.46** |
| 2 | **Chemistry: new shape or loaders** (fewer addresses, narrower words, loader-aware rounding) | The public moved +9.9% in 10 days, all from shape and loader changes, not weights. Loaders are 8,178 of the 14,840 Toffolis, i.e. 2 × 4,089, so 55%. A loader-bound argument says public loaders are near 2√(mL) for *random* words, and rounding chosen for the loader is untried (757030, author's claim). | ×1.10 to ×1.25 (unproven) | ×1.032 to ×1.077 | 26.18 to 27.32 |
| 3 | **DQI: Toffolis per shot** (cheaper uncompute of the in-place GF(2¹²) chains; shared products) | Uncompute is "roughly half of the mean Toffoli count" (755302, author-reported). The public found only ×1.011 in 10 days. Width is capped at 52 spare qubits. | ×1.05 to ×1.10 | ×1.016 to ×1.032 | 25.78 to 26.18 |
| 4 | Swap in today's best public parts (`9c708f7b`, `eae396ba`, `f21b8a87` npz plus V22 loaders) | Observed LB 25.40482 | ×1.0046 (product) | ×1.0015 | 25.405 |
| — | Closed: chemistry λ_eff refits; DQI width below 840; exploiting unreachable invsqrt x (about 0.08%, [753451](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753451)); R-reset of input-dependent scratch, which the verifier already catches (755302) | | | ≤ ×1.001 | |

For scale, gold today (rank 15) is 26.80, which needs ×1.057 on our blend. Lever 1 at mikelou1's level alone gets
26.46. Lever 1 plus a ×1.10 chemistry gain gets 27.31.

## Relevant primary sources

- **[arXiv:2610.23832](https://arxiv.org/abs/2610.23832)**, the task paper (Khattar, Rubin, Gidney, Broughton, Low,
  Babbush, Oct 2026). Verifiable structure means reversible logic, phase operations and measurement-based
  uncomputation, certified by fuzz testing. It is the design brief for all three tasks.
- **[arXiv:2510.10967](https://arxiv.org/abs/2510.10967)**, *Verifiable Quantum Advantage via Optimized DQI Circuits*
  (Khattar, Shutty, Gidney et al.). The Extended Euclidean Algorithm is the bottleneck. It uses implicit Bézout
  access and in-place designs, with leading-order space 2nb. This bears on lever 3.
- **[arXiv:2407.17966](https://arxiv.org/abs/2407.17966)**, *Rise of conditionally clean ancillae* (Khattar &
  Gidney, Quantum 2025). Already worth −294 Toffolis per shot in public DQI. The same idea applies to the invsqrt
  shifter and the chemistry loaders.
- **[arXiv:2611.09822](https://arxiv.org/abs/2611.09822)**, *Discovering formally verified quantum programs with AI*
  (DeepMind and Google, Nov 2026). Lean-verified agent loops. The abstract reports 11× on "simulations relevant to
  fusion" against its own baselines, with no mapping to this competition's references. Read the body before relying
  on it for invsqrt.
- **[arXiv:2308.12352](https://arxiv.org/abs/2308.12352)**, Rubin et al. The stopping-power source of the invsqrt
  task.

## What is not established

- **Any top team's per-task ratio.** The split of the ×1.57 between tasks is inference by elimination, not data.
- **mikelou1's 25.1× invsqrt.** It is author-reported only: no file, no method, not verified.
- **The chemistry loader bound.** The "~2√(mL)" figure and "public loaders are near it" are a forum author's claims
  for random words. Whether co-designed words beat it is open.
- **DQI's "uncompute ≈ half of Toffolis".** The figure is author-reported, and not checked against the `eae396ba`
  profile.
- **Eligibility of the public chemistry frontier.** Some accepted ANGLE packets are correct only on reachable states
  ([757522](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757522)), and
  there is no host ruling. maoxin123's V22 says it resets only provably-|0⟩ qubits, but we have not audited that.
- **Milestone-prize mechanics.** It is unsettled whether the scoring notebook itself must be public
  ([758095](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758095)).
- **arXiv:2603.28846.** It is linked from the competition pages but was not opened.

## Recommendation: which task to attack and why

**Attack `invsqrt` first. Keep chemistry second, and only through a new shape or loader co-design. Leave DQI last.**

1. **invsqrt has evidence of private headroom.**
   - It is the only task with a disclosed number above the public frontier: +13%, author-reported.
   - The public line moved +7.5% in 10 days, so it is not exhausted.
   - A qubit is worth 7× more there than in the other tasks.
   - The circuit is small (62 qubits, about 1.5k Toffoli-equivalents), so iteration is the fastest of the three.
   - Concrete open sub-problems with stated bounds: the shifter (high 40s against high 30s), the 62-qubit floor
     question, and MBU of the phase-stage uncompute.
2. **Chemistry gains need a structural change, not tuning.**
   - λ_eff is at the floor and the table fit is at the shape floor.
   - The public's +9.9% came from shape and loader changes.
   - Our prior loader-floor work (verdicts `r3`/`r4` in
     [`docs/experiments/verdicts.yaml`](../docs/experiments/verdicts.yaml)) says J=5 word-reading loaders cannot reach
     the band. So only a different table shape, or word rounding chosen for the loader, is worth the time.
   - It also carries an open eligibility ruling (757522).
3. **DQI is the most saturated task.** Its width is capped by the input register, and the public found ×1.011 in 10
   days.

**Before any of that:** swap in today's best public parts for +0.15% (25.40482 observed). This is a free submission
and not research.

**What changed against our plan**
([`docs/specs/2026-10-31-quantum-compilation-design.md`](../docs/specs/2026-10-31-quantum-compilation-design.md)):

- **The plan recommended `dfthcsos` for skill fit, and is now revised to `invsqrt` first.** Chemistry's continuous
  levers (λ, refits) are now proven closed. The only third-party evidence of headroom beyond public is on invsqrt.
- **The plan's DQI case ("888 → 800 is ~23% off the task") is wrong.** The verifier's 840-qubit syndrome register is
  the input, so the width can never go below 840, and the DQI width lever is capped at ×1.128 in theory.
- **The leaders' gap has held at ×1.16–1.17 since 2 Nov.** Assembling public parts will not close it, so our next
  points must come from our own compression.
