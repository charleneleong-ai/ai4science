# ai4science

AI-for-Science applications and research notes, organised around one through-line:

> **The generator is the commoditised part; the verifier is where reliability — and the value — actually lives.**

A fast generative/predictive model is increasingly a commodity. What converts a *plausible* output into a *trustworthy* one is an independent verifier — a deterministic check, a calibrated uncertainty band, a ground-truth oracle — that knows when the model can be trusted and when to defer. Across five projects this same "verification-first" pattern recurs:

| Project | Generator (imagination) | Verifier (the asset) |
| --- | --- | --- |
| [`nemotron-reasoning-challenge`](../nemotron-reasoning-challenge) | reasoning policy (Nemotron 30B) | RLVR — machine-checkable boxed answers (GSPO) |
| [`laguna-finetune`](../laguna-finetune) | agentic coding policy (Laguna XS.2) | reward-model-free RL — tests pass / render-diff / exact-match |
| [`orak-hackathon` / TGAER](../orak-hackathon) | symbolic reasoning agent | game env / world model as verifier |
| [`robotics_world_models`](../robotics_world_models) | learned world model (TD-MPC2, DreamerV3) | classical planner / when-to-trust crossover |
| [`smart_city_foundation_model`](../smart_city_foundation_model) | coupled weather→energy→demand forecast | conformal calibration — guaranteed uncertainty band |

This repo collects the applied work and research grounding where that thesis met two 2026 opportunities:

- [`materialhack/`](materialhack) — **MaterialHack** (Nucleate UK × ARIA, London, Jun 26–28 2026). AI track: protein–metal binders for critical-mineral recovery. The verifier thesis applied to de-novo protein design under extreme (out-of-distribution) conditions — built as [`touchstone/`](touchstone), a generator-agnostic metal-binder verifier.
- **Orbital Industries / Orb** — the atomistic-simulation foundation-model direction: the verifier thesis applied to materials, a trust layer on a fast learned surrogate for quantum reality. Now folded into touchstone as the [`MLIP UQ trust-layer`](touchstone/docs/specs/2026-06-30-mlip-uq-trust-layer.md) design spec, with a first step shipped (relaxation-health gate + OrbMol backbone).

The one-line pitch that ties it all together:

> A verifiable world model for the planet — weather, energy, and movement turned into trusted decisions. The same architecture powers embodied reasoning and materials foundation models: a fast generative world model, gated by a verifier that knows when to trust it.

---

## For my own learnings — an intro to AI4Science modelling

A long-living learning log. I come from the RL / verifier / world-model side, not classical computational chemistry or biology, so these are my notes-to-self on how the field is shaped and where my existing skills transfer. I'll keep adding to it.

This section is the evergreen overview; **dated, distilled entries live in [`learning-log.md`](learning-log.md)** (maintained by the `ai4science-log` skill — ground-first, then distil what transfers).

### The one pattern under all of it
AI4Science keeps doing the same move: **replace an expensive physical ground truth with a fast learned surrogate, then use it in a loop.**

| Domain | Expensive ground truth | Fast learned surrogate |
| --- | --- | --- |
| Materials / atoms | DFT (quantum, ~O(N³)) | MLIP foundation models (Orb, MACE, CHGNet, MatterSim) |
| Molecules / chemistry | high-level QM | OrbMol, long-range-electrostatics potentials |
| Proteins / bio | wet-lab assay, crystallography | structure/co-fold (AlphaFold3, Boltz-2), design (RFdiffusion, ProteinMPNN) |
| Embodied / robotics | real-world rollout | world models (TD-MPC2, DreamerV3, World-Action Models) |
| Earth / geo | sensor reality | weather FMs (Aurora), time-series FMs (Chronos-2) |

The catch — and my way in: **a surrogate is fast but goes out-of-distribution and is then confidently wrong.** So the field's real bottleneck isn't a bigger generator, it's a *trust layer* — uncertainty quantification, OOD detection, calibration, active learning. That's the verification-first thesis above, and it's domain-agnostic.

### The transferable techniques (what's worth getting fluent in)
- **Foundation-model pretraining** — self-supervised recipes like denoising/diffusion (Orb's atom recipe; same idea as generative models).
- **Fine-tuning / transfer** — LoRA, multihead replay; "initialisation/hparams > method" ([2606.12704](https://hf.co/papers/2606.12704)). LoRA on Nemotron/Laguna/Chronos is directly portable.
- **Uncertainty & OOD** — deep ensembles, Bayesian-by-disagreement ([2510.03046](https://hf.co/papers/2510.03046)), energy-based OOD ([2010.03759](https://hf.co/papers/2010.03759)), conformal calibration (used in the city twin). Caveat: naïve UQ ≠ reliable OOD ([2012.05329](https://hf.co/papers/2012.05329)).
- **Active learning / autoresearch** — spend expensive ground truth only where the surrogate is uncertain. This *is* sample-efficient autoresearch (cf. my `autoresearch.py`).
- **RL with verifiable rewards** — RLVR/GSPO; the reward *is* the verifier (my Nemotron/Laguna work).

### Domain vocab I'm still building (gaps to close)
- The actual physics: DFT, exchange-correlation functionals, potential-energy surfaces, **conservative vs non-conservative forces** and why energy conservation matters for MD stability.
- **Equivariant architectures** (MACE, e3nn, the equivariance/conservatism/sparsity trade-off) — Orb's bet that you can lean away from hard equivariance.
- Hands-on protein design tooling (RFdiffusionAA, ProteinMPNN, co-folding) and metal coordination chemistry.
- Benchmarks as ground for claims: **Matbench Discovery** (crystal stability), and the benchmark-vs-experiment gap ([UniFFBench, 2508.05762](https://hf.co/papers/2508.05762)).

### Reading list (grounded June 2026)
- **Orb** [2410.22570](https://arxiv.org/abs/2410.22570) · **Orb-v3** [2504.06231](https://arxiv.org/abs/2504.06231) — the materials FM I'm centring on.
- **When to Trust Imagination** [2605.06222](https://hf.co/papers/2605.06222) — the verifier thesis, cleanest single instance.
- **World Model for Robot Learning: a survey** [2605.00080](https://hf.co/papers/2605.00080) — the embodied-substrate map.
- **Towards Agentic Intelligence for Materials Science** [2602.00169](https://hf.co/papers/2602.00169) — autoresearch/agent framing.

> Working principle for everything in here: *establish the verifier before trusting the generator.* Same discipline as a baseline before an ML experiment — know what "right" means before you optimise toward it.

---

Grounded against research current as of **June 2026**.
