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

This repo collects the application materials and research grounding where that thesis met two 2026 opportunities:

- [`orbital/`](orbital) — **Orbital Industries** (the Orb / atomistic-simulation foundation-model company). ML Engineer application. The verifier thesis applied to materials: a trust layer on a fast learned surrogate for quantum reality.
- [`materialhack/`](materialhack) — **MaterialHack** (Nucleate UK × ARIA, London, Jun 26–28 2026). AI track: protein–metal binders for critical-mineral recovery. The verifier thesis applied to de-novo protein design under extreme (out-of-distribution) conditions.

The one-line pitch that ties it all together:

> A verifiable world model for the planet — weather, energy, and movement turned into trusted decisions. The same architecture powers embodied reasoning and materials foundation models: a fast generative world model, gated by a verifier that knows when to trust it.

Grounded against research current as of **June 2026**.
