# The DF-THC one-norm is 94% a constant: why no factorization refit can lift this score

**Date:** 2026-11-02 · **For:** the `dfthcsos` task, after gate G1 cleared chemistry for
single-seed comparison · **Status of our entry:** blend 23.18119, chemistry cost 337,780,865
at Λ = 299.82, λ_eff = 48.293.

**Headline:** the literature's strongest one-norm lever (joint symmetry-shift optimisation) is
frozen by the verifier, and the one it leaves open — an L1-regularised refit of the THC factors
— is worth at most **+0.243 blend** against the **+1.82** needed for 25.0, because `lambda_eff`
is **93.9% a fixed constant of the reference Hamiltonian**. Proven below, from the verifier's
own source, and verified numerically to machine precision. The refit is **DISCARDed before
launch**; the chemistry score lives in `T` and `Q`.

## Why this was asked

The chemistry cost is `lambda_eff * T * Q` with `lambda_eff = sqrt(E_gap * (2Λ - E_gap))`, and
measurement (below, and in `docs/experiments/`) showed the submitted factorization pinned
against every cheap knob: both Frobenius residuals at 99.999% of the 20 mHa cap, `E_gap` 1e-6
above its strict floor, `brot` already down at 8, and Λ 95.7% `lambda_h2`. Nothing was gettable
by tuning, so the question became what the field's best *fitting* methods are and which of them
this verifier permits — and, as it turned out, whether a better fit could matter at all.

One earlier reading is corrected below: a `W` scan that appeared to show `error_1` holding `W`
at a sharp local minimum had held the submitted one-body rows fixed — see *Where the chemistry
score actually lives*.

## The method the task is built from

[Lee, Berry, Gidney, Huggins, McClean, Wiebe and Babbush, *Even more efficient quantum
computations of chemistry through tensor hypercontraction*, PRX Quantum 2 (2021)
(arXiv:2011.03494)](https://arxiv.org/abs/2011.03494) is the parent method: THC-factorize the
Coulomb operator, transform it to an isospectral diagonal form in the non-orthogonal basis the
THC factors define, and qubitize that. Toffoli complexity is `O(N*lambda/eps)`, which is why
the one-norm, not the rank, is what the score is most sensitive to — exactly the structure the
verifier reproduces.

## Where the field gets its one-norm reductions

Two distinct levers, and the distinction decides this experiment.

**Lever 1 — the symmetry shift.** Subtract a function of the electron-number operator,
`H~ = H + B(N_e - n_e)`. Since `N_e` commutes with `H`, the target sector's spectrum is
untouched and the energy is recovered by an a-posteriori correction, but the one-norm falls.
[Loaiza and Izmaylov's Block-Invariant Symmetry Shift (BLISS, arXiv:2304.13772)](https://arxiv.org/abs/2304.13772)
generalises it to an operator that need not be an exact symmetry;
[Patel et al. (arXiv:2409.18277)](https://arxiv.org/abs/2409.18277) find its global optimum by
linear programming, noting that earlier nonlinear optimisation could miss it.
[Rocca et al., *Symmetry-compressed double factorization* (arXiv:2403.03502, JCTC 20, 2024)](https://arxiv.org/abs/2403.03502)
optimise the shift prefactors *jointly with the factorization*, and report that zeroing those
prefactors at `N_DF = 5N` sends the one-norm back up to 145.8 Ha.
[Follow-up work (arXiv:2412.01338)](https://arxiv.org/abs/2412.01338) claims a further 25%
reduction in the block-encoding constant by optimising shifts and factorization together, and
[BLISS-THC (arXiv:2501.06165)](https://arxiv.org/abs/2501.06165) carries the shift into THC
specifically — so the shift is not a side technique but the axis the recent work is built on.

**This lever is closed to us.** The verifier's `load_inputs` calls `load_fixed_shifts()` and
*overwrites* whatever `H1_I` and `Hsym_lower` a submission carries, before `S` and `beta_1` are
derived — the shift is always the baseline's. The competition has frozen the single technique
the recent literature leans on hardest. Any plan that assumes a BLISS-style gain is planning
for a rejected submission.

**Lever 2 — the fit itself, at fixed shift.** This is the one the verifier leaves open,
and the rest of this report prices it. The established recipe,
from [Goings et al. on cytochrome P450 (PNAS 119, 2022)](https://www.pnas.org/doi/10.1073/pnas.2203533119),
is: initialise from a symmetric canonical-polyadic decomposition of the Cholesky vectors, then
run **L1-regularized optimization of the THC factors**. The regularization exponent is the
choice — 1 for L1, 2 for L2 — and the symmetry-compressed DF paper notes prior L2 work already
reduced λ relative to plain XDF and THC. BLISS-THC's contribution on this axis is to optimise
the ℓ2 fit error *and* the 1-norm together, trading rank against λ, which is precisely our
constrained problem with the shift term deleted.
[*Benchmarking least-squares THC techniques* (arXiv:2608.17885)](https://arxiv.org/abs/2608.17885)
is the most recent survey of the fitting step itself — the route that would have been
followed, had it been worth following.

## The objective, and why it is already exhausted

The objective is written for us by the verifier's own source:

```
minimise   lambda_h2(w, W) = 1/4 * sum_(rho,c) ( ||w^(rho c)||_1 + |W^(rho c) + sum_b w_b^(rho c)| )^2
subject to error_2 = || sum_(rho c) X (x) X - shift - tei ||_F <= 20 mHa,
           error_1 <= 20 mHa,  E_gap > 3.914898,  Lambda > 49.0
where     X^(rho c) = sum_b w_b^(rho c) u_b^(rho) u_b^(rho)^T
```

That is an **L1²-penalised THC fit with the symmetry shift held fixed** — lever 2 exactly. The
plan was to run it: gradient descent on `(u, w)` from the incumbent, constraints as penalties,
`brot = 8` by straight-through rounding.

**That experiment is closed before launch.** Writing the score's first factor out in full shows
the reachable range is 3.1%, and the proof is three lines of algebra on the verifier's own
quantities.

`lambda_eff = sqrt(E_gap · D)` with `D = 2·Lambda − E_gap`, `Lambda = P + M + lambda_h2` and
`E_gap = c0 + 2M + ½ΣW²`, where `P` and `M` are the positive and negative eigenvalue mass of
`h1_prime`. Substituting, **`M` cancels**:

```
D = 2P + 2·lambda_h2 − ½ΣW² − c0  =  2P + 2Δ − c0,     Δ := lambda_h2 − ¼ΣW²
```

`Δ ≥ 0` identically, not incidentally: termwise `A + B ≥ |W|` where `A = ‖w^ρc‖₁ ≥ 0` and
`B = |W^ρc + Σ_b w_b^ρc| ≥ |W^ρc| − A`, so every `(A+B)² − W²` term is non-negative. Measured
on the incumbent: all 480 terms positive, smallest `4.1e-5`, and `Δ` recovered to `2e-14`
([`TestTheFloorOnTheScore`](../tests/test_chemistry.py)).

So `D ≥ −c0`, and `c0 = E_upper − ecore + n·β₁` is **fixed**: `UPPER_ENERGY = −22140.27335`
and `ecore` are constants of the shipped verifier's `scoreinfra/reference.py`, and `β₁` comes
from its `load_fixed_shifts()`. It evaluates to **−559.132473**
([`FEMOCO_C0`](../src/quantum_compilation/cost.py)), which is **93.86% of the incumbent's `D` of
595.73**. With `E_GAP_FLOOR = 3.914898` likewise hardcoded
([`min_lcu_lambda`](../src/quantum_compilation/cost.py)):

| | value |
| --- | --- |
| `lambda_eff` floor = `sqrt(3.914898 · 559.132473)` | **46.786180** |
| incumbent `lambda_eff` | 48.293182 |
| best conceivable cut | **3.121%** |
| blend gain, at +0.0778 per 1% of a task's cost | **+0.243** |

And 46.786 is unreachable, not merely a bound: it needs `P = 0` *and* `Δ = 0`, and `Δ = 0`
forces `w ≡ 0`, which destroys `error_2` entirely. Halving `Δ` — already an optimistic ask
against a fit whose author reports a weights-only reoptimisation bought 0.005% — yields
`lambda_eff` 47.60, **+0.11 blend**.

A perfect refit therefore delivers 13% of the **+1.82** the blend needs, and a realistic one
about 6%.

**Verdict: DISCARD, before any compute was spent.**

## Where the chemistry score actually lives

`cost = lambda_eff · T · Q = 48.293 · 16,266 · 430 = 3.3778e8`, against a reference of 3.3556e9
for a ratio of **9.93** — the weakest of our three tasks (`invsqrt` 20.62, `dqi` 60.80), so the
headroom is here if it is anywhere. But it is not in `lambda_eff`. With that factor floored
within 3.1%, the chemistry cost is **`T · Q` and essentially nothing else**, and both are
circuit-compilation quantities — the Toffoli count of the ANGLE loaders and the 430-qubit
footprint — not properties of the factorization.

This also reframes the two saturated error budgets. Both residuals sit at 99.999% of the 20 mHa
cap, which is not a fit straining to be feasible: `error_1` compares the submitted one-body rows
against the exact eigenvectors of `h1_prime` and is **decoupled from every scored quantity**
(the eigenvalues come from the exact `h1_prime`, never from the submitted rows). Spending that
budget buys coarser rotations, which buys Toffolis. The incumbent has already sold both error
budgets for `T`. That is the market to trade in, and the next measurement must price it before
committing: **what does one Toffoli of the ANGLE loaders cost in error budget, and how many are
left to buy?**

The 20.4% that `25.0` needs off the product of the three costs is not in the chemistry
factorization. Nothing in this sweep says it is not in the loaders.

## Not established

- **No FeMoCo λ value from the literature is quoted here**, because none of the sources found
  states one in a form comparable to this verifier's `Lambda`, which includes the fixed
  symmetry shift and this task's own normalisation. The Lee et al. FeMoCo figures that did
  surface (≈5.3×10⁹ Toffolis, 2,142 logical qubits for a 108-qubit representation) are for a
  different resource model than the one scored here and must not be compared to our 16,266.
  Establishing a comparable published floor for Λ would need the primary papers read in full.
- **Whether `T` and `Q` have any headroom at all is unmeasured.** The `lambda_eff` floor says
  where the chemistry score is *not*; it does not say the loaders can be improved. Pricing one
  Toffoli in error budget is the next measurement, and it may come back at zero.
- The 21% BLISS figure is quoted from a survey, not reproduced, and is moot here anyway.
- The floor proof assumes the shipped verifier is the scored one. It is the dev-set verifier,
  and the three constants it rests on (`UPPER_ENERGY`, `ecore`, `E_GAP_FLOOR`) are hardcoded
  there. If the leaderboard scorer differs, the floor moves with `c0`.
