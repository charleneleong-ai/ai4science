# AI4Science learning log

Dated, distilled entries — newest first. Each entry: what it is / why it matters / how it transfers to my work / what to learn next. Maintained by the `ai4science-log` Claude skill (or by hand, same format). Evergreen overview + reading list live in [`README.md`](README.md#for-my-own-learnings--an-intro-to-ai4science-modelling).

<!-- NEW ENTRIES GO DIRECTLY BELOW THIS LINE (newest first) -->

## 2026-06-16 — Orb (universal interatomic potential)

**Source:** [Orb v1 — arXiv 2410.22570](https://arxiv.org/abs/2410.22570) · [Orb-v3 — arXiv 2504.06231](https://arxiv.org/abs/2504.06231) (grounded 2026-06-16)
**What it is:** A foundation model for atoms — input a 3D structure, output energy + forces; a fast drop-in for DFT in geometry optimisation, MD, Monte Carlo. Recipe: diffusion/denoising pretraining → supervised NNP. v3 charts an equivariance/conservatism/sparsity Pareto frontier and reaches the mesoscale (>10k atoms).
**Why it matters:** It's the "fast learned surrogate for quantum reality" — the verifier/world-model layer of the materials stack, and Orbital's crown jewel.
**How it transfers:** Same surrogate-in-a-loop pattern as my world-model and forecasting work. My entry point is the trust layer (UQ/OOD/conformal) that decides when to trust Orb vs. recompute DFT. Learn: conservative vs non-conservative forces, why energy conservation matters for MD stability.
**To learn next:** equivariant architectures (MACE/e3nn); how `equigrad` works.

## 2026-06-15 — When to Trust Imagination (World-Action Models)

**Source:** [arXiv 2605.06222](https://hf.co/papers/2605.06222) (grounded 2026-06-15)
**What it is:** World-Action Models imagine future frames + actions but execute a fixed chunk blindly. Reframes execution as a *future–reality verification problem*; adds FFDC, a lightweight verifier that compares predicted vs. real observation each step and replans when they diverge.
**Why it matters:** They didn't improve the world model — they bolted a *verifier* on the same imagination, and real-world success went 45% → 80%. Cleanest single instance of "the verifier is the asset."
**How it transfers:** Direct analog of RLVR (reward = verifier) and conformal forecasting (band = verifier). The materials version is a trust layer on Orb; the harder twist there is that ground truth (DFT) is expensive, unlike a robot's free camera frame.
**To learn next:** how FFDC trains its verifier on synthetic negatives — does that idea port to UQ for surrogates?
