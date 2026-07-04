# MLIP uncertainty / OOD trust layer — design

**Date:** 2026-06-30
**Status:** research grounding (future direction, feeds phase 3)
**Builds on:** [`2026-06-26-mlip-and-structural-verifier-design.md`](2026-06-26-mlip-and-structural-verifier-design.md),
[`2026-06-28-rlvr-boltzgen.md`](2026-06-28-rlvr-boltzgen.md)

## Why

The MLIP tier today relaxes the metal cluster and ranks ΔE, but it has **no notion
of when to trust the MLIP itself**. A foundation interatomic potential is a fast
learned surrogate for quantum reality; in a relaxation or an MD rollout it commits
to that surrogate and, when the configuration wanders out-of-distribution, it can be
*confidently wrong*. That is the same blind spot touchstone exists to catch in
generators — only now the unreliable model is one of touchstone's own tiers.

The unlock is the same pattern touchstone already ships everywhere else: **the
generator is the commodity, the verifier is where reliability lives.** Applied to
the MLIP, that means a calibrated **uncertainty / OOD signal on the potential's own
prediction** that decides where the MLIP's energy can be believed and where an
expensive higher-fidelity check (DFT) must be paid for.

## Origin — the MaterialHack brief (the OOD thesis)

This trust layer is the direct descendant of the project's founding framing (MaterialHack
2026, Nucleate UK / ARIA, AI track — *"design protein–metal interactions for critical-mineral
recovery"*). Of the track's three questions, **Q2 — "accurately predict protein behaviour in
extreme environments" — is the verification question in disguise.** Real critical-mineral recovery
(bioleaching ores, e-waste, brines, mine tailings) happens in **acidic (pH 1–2), hot, high-salinity
leachate** — exactly where a model trained on benign-condition data goes out-of-distribution and
**confidently wrong**. The contrarian move — versus teams that just maximise predicted binding
affinity — was to **design for robustness and know when you can't be trusted.**

Two commitments from that brief carry straight into touchstone and this layer:
- **The verifier must be independent of the generator**, and can be validated *in silico* — via a
  physics oracle (co-fold, MLIP) and a held-out slice of the existing experimental record
  (MetalPDB / BioLiP), no wet lab required. This MLIP-UQ layer is that principle turned on
  touchstone's *own* MLIP tier.
- **Extreme-condition inputs need an explicit OOD flag.** touchstone's `--stress` map
  (`neutral` / `leachate` / `low_pH`) already probes the operating envelope at the geometry tier;
  this layer extends "flag when you can't be trusted" to the physics tier. (Founding brief archived
  here; the standalone `materialhack/` note was folded in and removed.)

## The structural mapping

The robotics result this borrows from is *When to Trust Imagination: Adaptive Action
Execution for World-Action Models* ([arXiv 2605.06222](https://arxiv.org/abs/2605.06222)):
a world model imagines future frames, executes a fixed chunk blindly, and a
lightweight verifier (FFDC) compares predicted-vs-real each step and replans when
reality diverges. They did not improve the world model — they bolted a verifier on
the same imagination and real-world success went 45% → 80%. The intelligence that
mattered was in *checking* the future, not generating it.

| *When to Trust Imagination* (robotics) | touchstone MLIP tier |
| --- | --- |
| World model imagines future frames + actions | MLIP imagines energy + forces (the PES) |
| Blindly executes a fixed action chunk | Blindly rolls out MD / ranks ΔE |
| "Reality deviates from imagination" | Trajectory drifts **out-of-distribution** |
| FFDC verifier → confidence eₜ ∈ [0,1] | **UQ / OOD score** on the MLIP prediction |
| eₜ < 0.5 → replan, re-ground in real obs | Score high → **trigger a DFT recompute** |
| Adaptive chunk size from consistency | Spend DFT **only where the MLIP is unreliable** |
| Verifier trained on synthetic negatives | UQ calibrated on OOD / high-error configs |

Both are *adaptive compute allocation gated by a verifier*: cheap surrogate where
trustworthy, expensive ground truth only where it breaks.

**The sharp disanalogy.** For a robot, reality is cheap — the camera gives a free
ground-truth frame every step. For an MLIP, ground truth (DFT) is the *expensive*
thing, so the verifier must estimate its own unreliability **without observing
reality** and only then pay for the check. The UQ head plays the role FFDC got for
free from the camera. This is the harder version of the problem, and it is exactly
what makes it the right verifier-first investment for touchstone.

**The hook already exists.** Touchstone is not missing a trust layer — every
verdict carries an `ood` flag ([`core.py` `Verdict.ood`](../../src/touchstone/core.py)),
the geometry tier ships a real OOD probe ([`geometry/ood.py`](../../src/touchstone/geometry/ood.py),
extreme-leachate perturbation), and the MLIP tier already defers on `site_drift >
ood_drift`. What this layer adds is making the MLIP's *own* uncertainty drive that
defer, rather than a geometric drift heuristic standing in for it.

## Why this sharpens the existing tiers

- **MLIP binding tier** — the ΔE ranking is explicitly ranking-only today because
  absolute solution-phase binding energies from MLIPs are unreliable. An OOD score
  turns "unreliable everywhere, so rank-only" into "reliable here, suspect there,"
  letting borderline designs escalate to DFT instead of being silently mis-ranked.
- **MLIP-MD (phase 2 structural axis)** — long MD needs **energy conservation**, so
  the backbone must use **conservative forces** (F = −∂E/∂x, gradient of one scalar
  energy) rather than a direct-force head, which is faster/sparser but carries no
  conservation guarantee and drifts over ns-scale rollouts. A UQ signal additionally
  flags when the trajectory has left the training manifold mid-rollout.
- **RLVR-for-BoltzGen autoresearch loop** — a verifier-gated reward already drives
  reject-sampling fine-tuning. A calibrated UQ head makes that loop **sample-
  efficient**: spend DFT (the costly ground-truth reward) only on candidates where
  the MLIP reward is untrustworthy, not uniformly.

## Candidate mechanisms, cheapest first

- **Relaxation health — implemented.** The single relaxation the MLIP tier already
  runs exposes its own "did I settle?" signal for free: whether LBFGS reached `fmax`
  within the step budget, and the final `|F|max`. A non-converged relaxation has read
  its drift/CN off a non-equilibrium geometry, so the site can't be vouched for →
  `defer`. No second model, no new dependency, no calibration data needed — the
  materials analog of FFDC's per-step confidence. Shipped on
  [`MLIPVerifier.verify`](../../src/touchstone/physics/mlip.py) /
  [`SiteRelaxation`](../../src/touchstone/physics/mlip.py).
- **Ensemble disagreement** ([Bayesian E(3) potential,
  arXiv 2510.03046](https://arxiv.org/abs/2510.03046)) — the most defensible UQ
  signal: run a second backbone on the same cluster and `defer` when the two disagree
  on the CN-call / ΔE / geometry. The pairing matters on two axes — deployment cost and
  how *differently* the two models fail:
  - **MACE ↔ Orb is the preferred in-process pair.** [orb-models](https://github.com/orbital-materials/orb-models)
    declares **no `e3nn`** (Orb is non-equivariant), so it co-resides with MACE's
    `e3nn==0.4.4` in one env — verified by resolving both together: `mace-torch 0.3.16`
    + `orb-models 0.7.0` + `torch 2.12.1` + `e3nn 0.4.4`, clean, no conflict. Apache-2.0,
    and the *most architecturally diverse* pairing (equivariant MP-crystal vs
    non-equivariant OMol25-molecular), so they fail differently — what a disagreement
    signal needs. [OrbMol](https://hf.co/Orbital-Materials/OrbMol) is also better-posed
    for this target: it takes total **charge + spin** (the cluster is a charged
    open-shell Ni²⁺ site), trains on OMol25's transition-metal complexes, and v2 adds
    long-range electrostatics — exactly the regime where short-range MACE's ΔE is weak.
    Unproven in touchstone, though: Orb must first reproduce the GFN2/MACE/UMA Ni-cluster
    CN call (the benchmark UMA passed) before it's trusted as the second opinion.
  - **UMA stays the separate-env third opinion.** It already earned its keep
    corroborating the CN 4→3 labile-donor call, but `fairchem` pins `e3nn 0.6` against
    MACE's `0.4.4`, so it needs its own env (subprocess orchestration) and is gated
    (FAIR AUP). Keep it as the heavier check to spin up on demand or to break a
    MACE↔Orb tie — not dropped, just no longer the default second backbone.
- **Energy-based OOD** ([arXiv 2010.03759](https://arxiv.org/abs/2010.03759)) —
  appealing because the MLIP already emits an energy, but a raw potential energy is
  not a calibrated OOD score the way a classifier's logits are; it needs a
  per-atom-energy or force-residual outlier model and a reference distribution. A
  research step, not a first-pass freebie.

**Counter-evidence to respect (so the calibration is honest):**
- naïve UQ ≠ reliable OOD detection — *Know Your Limits*
  ([arXiv 2012.05329](https://arxiv.org/abs/2012.05329)).
- benchmark-vs-experiment gaps are real for MLIP force fields — UniFFBench
  ([arXiv 2508.05762](https://arxiv.org/abs/2508.05762)). A UQ score calibrated only
  against DFT can still be wrong about the lab.

## Sequencing

Folds into **phase 3 (calibration)** of the MLIP design.

1. **Relaxation-health gate — done.** Non-convergence / high `|F|max` → `defer` on the
   MLIP verdict. Single backbone, zero new deps.
2. **Ensemble disagreement** — the `orbmol` backbone is wired (charge/spin plumbed
   through to `atoms.info`, `[orb]` extra, fail-fast when charge/spin are missing);
   remaining is the MACE ↔ Orb disagreement comparison itself, gated on Orb first
   passing the Ni-cluster CN benchmark UMA already passed. UMA stays the separate-env
   third opinion.
3. **Calibrated UQ head** against DFT / high-error configs — the research-program
   version, gated on wet-lab outcome data touchstone does not yet have.
