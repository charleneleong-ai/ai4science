# Cutting the DF-THC one-norm: what the literature does, and which half of it this competition forbids

**Date:** 2026-11-02 · **For:** the `dfthcsos` task, after gate G1 cleared chemistry for
single-seed comparison · **Status of our entry:** blend 23.18119, chemistry cost 337,780,865
at Λ = 299.82, λ_eff = 48.293.

## Why this sweep was run

The chemistry cost is `lambda_eff * T * Q` with `lambda_eff = sqrt(E_gap * (2Λ - E_gap))`, and
measurement (below, and in `docs/experiments/`) shows the submitted factorization is pinned
against every cheap knob: both Frobenius residuals sit at 99.999% of the 20 mHa cap, `E_gap`
sits 1e-6 above its strict floor, `W` is held by `error_1` at a sharp local minimum, and `brot`
is already down at 8. Λ is 95.7% `lambda_h2`, the two-body LCU one-norm. Nothing is gettable by
tuning; the only route left is a better *fit*. This asks what the field's best fitting methods
are, and which of them this verifier permits.

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
specifically. A 2025 block-encoding survey puts BLISS's effect on THC LCUs at **roughly 21% on
average**.

**This lever is closed to us.** The verifier's `load_inputs` calls `load_fixed_shifts()` and
*overwrites* whatever `H1_I` and `Hsym_lower` a submission carries, before `S` and `beta_1` are
derived — the shift is always the baseline's. The competition has frozen the single technique
the recent literature leans on hardest. Any plan that assumes a BLISS-style gain is planning
for a rejected submission.

**Lever 2 — the fit itself, at fixed shift.** This is what remains. The established recipe,
from [Goings et al. on cytochrome P450 (PNAS 119, 2022)](https://www.pnas.org/doi/10.1073/pnas.2203533119),
is: initialise from a symmetric canonical-polyadic decomposition of the Cholesky vectors, then
run **L1-regularized optimization of the THC factors**. The regularization exponent is the
choice — 1 for L1, 2 for L2 — and the symmetry-compressed DF paper notes prior L2 work already
reduced λ relative to plain XDF and THC. BLISS-THC's contribution on this axis is to optimise
the ℓ2 fit error *and* the 1-norm together, trading rank against λ, which is precisely our
constrained problem with the shift term deleted.
[*Benchmarking least-squares THC techniques* (arXiv:2608.17885)](https://arxiv.org/abs/2608.17885)
is the most recent survey of the fitting step itself and is the first thing to read before
choosing an optimiser.

## What this means for the plan

The objective is already written for us by the verifier's own source:

```
minimise   lambda_h2(w, W) = 1/4 * sum_(rho,c) ( ||w^(rho c)||_1 + |W^(rho c) + sum_b w_b^(rho c)| )^2
subject to error_2 = || sum_(rho c) X (x) X - shift - tei ||_F <= 20 mHa,
           error_1 <= 20 mHa,  E_gap > 3.914898,  Lambda > 49.0
where     X^(rho c) = sum_b w_b^(rho c) u_b^(rho) u_b^(rho)^T
```

That is an **L1²-penalised THC fit with the symmetry shift held fixed** — lever 2 exactly, and
nothing else. Parameters are `u` (8×143×76), `w` (8×143×60) and `W` (8×60), about 155k values,
against a 76⁴ integral tensor; gradient-based optimisation is tractable on CPU.

**First thing to try:** L1-regularised refit of `(u, w)` by gradient descent from the current
best factorization as the initialisation, with the constraints as penalty terms and the
quantisation to `brot = 8` handled by straight-through rounding. Starting from the incumbent
rather than from a CP decomposition is the deviation from the published recipe, and it is
deliberate: the incumbent already satisfies the constraints, so the optimiser starts inside the
feasible set and the question becomes only whether λ_h2 can be walked downhill from there.

**What would make this fail, stated up front:** the incumbent may already be at a local
optimum of the same objective — its author's notes say a weights-only reoptimisation bought
0.005%. If gradient descent from the incumbent moves λ_h2 by less than a percent, the answer is
that a *different initialisation* (CP decomposition of the Cholesky vectors, per the published
recipe) or a different shape is needed, and that is a larger experiment.

## Not established

- **No FeMoCo λ value from the literature is quoted here**, because none of the sources found
  states one in a form comparable to this verifier's `Lambda`, which includes the fixed
  symmetry shift and this task's own normalisation. The Lee et al. FeMoCo figures that did
  surface (≈5.3×10⁹ Toffolis, 2,142 logical qubits for a 108-qubit representation) are for a
  different resource model than the one scored here and must not be compared to our 16,266.
  Establishing a comparable published floor for Λ would need the primary papers read in full.
- **No code was located that performs this fit**, only papers describing it. The P450 paper
  provides resource *estimation* as an OpenFermion submodule, which is not the fitting routine.
  Whether `openfermion`, `pyscf` or the BLISS-THC authors ship a usable optimiser is unchecked.
- The 21% BLISS figure is quoted from a survey, not reproduced, and is moot here anyway.
