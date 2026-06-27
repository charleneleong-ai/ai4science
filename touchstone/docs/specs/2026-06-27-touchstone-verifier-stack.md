# Touchstone verifier stack — architecture

**Date:** 2026-06-27
**Status:** living spec (consolidates the stack as built across PRs #2/#3/#5/#6/#7/#8)
**Related:** [`2026-06-24-touchstone-design.md`](2026-06-24-touchstone-design.md) (original thesis),
the MLIP + structural-verifier scoping doc (PR #5)

## Thesis

The generator is a commodity; the **verifier is the asset**. A design flows
generator → `BinderDesign` → verifier stack → `Verdict`. The verifier never sees
*how* a design was made — only the `CoordinationSite` — which is what makes it
generator-blind, and what makes it the durable, reusable piece as generators churn.

## Data flow

```
Generators  (RFdiffusionAA · BoltzGen · Mock)
     │  emit
     ▼
BinderDesign { sequence, CoordinationSite, generator, source }   ← the only coupling point
     │  consumed by (blind to generator)
     ▼
VERIFIER STACK  ───►  Verdict { score (↑=better, ~[0,1]), trust, ood, reason }
                                     └─► label: trust / weak / defer
     │
     ▼
Pipeline:  rank() · design_and_rank() · selectivity_profile()
```

## The verifier stack — 4 stages, ≥2 independent methods each

Defense-in-depth: each stage is checked by two methods that fail differently, so a
design must satisfy all of them. Two independent verifiers agreeing is the core
reliability signal (it is also what makes the stack a hard-to-hack RL reward).

| stage | method A | method B | question | PR |
| --- | --- | --- | --- | --- |
| **Geometry** | `GeometryVerifier` (z-score vs reference) | `BondValenceVerifier` (Σ bond-valence ≈ formal charge?) | is the coordination plausible? | #5 |
| **Co-fold** | Boltz-2 | Chai-1 + AllMetal3D (via `CofoldCrossCheck`) | does an *independent* predictor agree? | #6 |
| **Physics (statics)** | xtb GFN2 | `MLIPVerifier` (MACE-MP / UMA) | does the site hold under relaxation? | #5 |
| **Dynamics** | xtb cluster-MD | `MLIPDynamicsVerifier` | does it survive 300 K MD? | #5 |

`CofoldCrossCheck` is predictor-agnostic — Boltz-2, Chai-1, and AllMetal3D plug in
via a `provider` callback, no verifier change.

## Pluggable seams ("swap, don't rewrite")

- **Reference oracle** (`ReferenceDistribution`): `MockReference` → `PDBReference` →
  `CSDReference` — all share `_JsonReference`; the verifier is unchanged whichever
  is plugged in (PR #7). CSD's metal–organic priors complement the PDB's protein
  sites.
- **MLIP backbone**: MACE-MP-0 (default, MIT) / UMA (gated) / EMT (test) — lazy,
  injectable, swappable (PR #5). Backbone chosen by benchmark vs the xtb GFN2
  optimum: MACE-MP 0.056 Å, UMA 0.083 Å; both independently make the same CN call,
  so it is a real signal, not a model artefact.
- **Verdict semantics**: `Verdict.defer()` centralises the trust/weak/defer contract
  across all verifiers; every score is higher=better in ~[0,1] so `rank()` and
  `selectivity_profile()` compose over any verifier.

## Coverage and honest gaps

Mapped to the wet-lab validation triad:

- **Binding** ✅ fully — the whole stack.
- **Thermostability** 🟡 *site-level* (MD survival of the coordination), **not**
  whole-protein Tm / global unfolding.
- **Expression** ❌ not covered — the genuine missing stage (solubility / aggregation
  / fold-confidence). The next vertical to add.

Calibration against real wet-lab outcomes (turning "trust" into a measured hit-rate)
is unbuilt — the verdicts are physically-grounded but not yet experiment-calibrated.

## Live trace — real BoltzGen designs (2026-06-27)

Ran the stack end-to-end on four **real BoltzGen** Ni-motif designs (His/His/His-Cys
theozyme sites → N/N/N/S donors), with their backbones, on the A100 (`mlip` env,
MACE-MP-0 float64, CUDA). This is the first run against generator output rather than
hand-built fixtures.

| design | CN · donors | geometry | bond-valence | MLIP relax | reward |
| --- | --- | --- | --- | --- | --- |
| `ni_motif_00` | 3 · N/N/S | WEAK 0.265 | TRUST (BVS 1.71) | DEFER (diverged, drift 9.1 Å) | 0.0 |
| `ni_motif_01` | 3 · N/N/S | WEAK 0.289 | DEFER (BVS 1.09, Δ0.91) | DEFER (RuntimeError) | 0.0 |
| `ni_motif_02` | 4 · N/N/N/S | TRUST 0.268 | **TRUST (BVS 2.00)** | DEFER (diverged, drift 7.9 Å) | 0.0 |
| `ni_motif_03` | 4 · N/N/O/S | TRUST 0.219 | TRUST (BVS 1.88) | DEFER (RuntimeError) | 0.0 |

**What works:** geometry + bond-valence discriminate cleanly and agree — the two
CN-4 designs (`ni_motif_02` His₃Cys, BVS 2.00 exact; `ni_motif_03`, BVS 1.88) read
as the real binders; the CN-3 pair reads weak/under-coordinated (BVS 1.71 / 1.09).
This is the stack doing its job on generator output without seeing the generator.

**What broke — the MLIP statics/dynamics tier is non-functional on backbone-bearing
designs as currently implemented.** `MLIPVerifier._cluster` extracts *metal +
donor atoms only*, dropping the protein backbone that physically holds the donors in
place. Relaxed as a free cluster, the donors disperse (drift 8–9 Å) and the
interaction energy blows up to ~10²⁷ eV, or the optimiser RuntimeErrors outright; the
200-step float64 MLIP-**MD** then hangs on the diverged/inf geometry, pinning all
80 GB of GPU at 0 % util (had to be killed). So every MLIP verdict DEFERs, which —
because consensus is min-like — drags **every** reward to 0.0 even for the BVS-2.00
design. The physics tier is actively *harming* the signal, not adding to it.

**Fix (next):** the cluster must carry the backbone — either include backbone atoms
and **position-restrain** them during relaxation/MD (so only the side-chain donors
relax against a fixed scaffold), or cap the residues and freeze Cα. Until then the
statics/dynamics tier should be **gated out** of the reward (geometry + bond-valence
+ co-fold only), not silently DEFERed into a 0.0. Tracked as the top physics-tier
TODO; the hang also argues for an MD pre-check that aborts on non-finite forces
instead of grinding for minutes.

## Forward — post-training (RLVR)

The verifier stack is, by construction, a **verifiable reward**: it emits a scalar
(`score` / `trust`) grounded in physics + experimental references, not in the
generator's own confidence. That makes it the reward signal for *post-training* a
generator — the `biometal_rlvr` direction.

```
Generator (BoltzGen / LigandMPNN)  ──generate──►  designs
        ▲                                              │
        │  policy update (RL / DPO / best-of-N)        ▼
        └──────────── reward ◄──── touchstone verifier stack (CSD-grounded)
```

### Which generator to post-train

BoltzGen is the natural target *here* — open-weight, all-atom, metal/ligand-aware,
the AI-track sponsor model, and it carries a **diversity edge** (higher Vendi score
than RFdiffusion / RFdiffusionAA, conditioned on target —
[Tamarind](https://www.tamarind.bio/blog/boltzgen-validated-de-novo-binder-design-for-diverse-targets),
[SeedProteo benchmark](https://arxiv.org/html/2512.24192)). Honest caveat: it is one
of a **cluster** (RFdiffusion(AA/3), BindCraft, PXDesign, AlphaProteo, SeedProteo),
and RFdiffusion/BindCraft are the *most-validated* — so "best candidate for
post-training" is about open weights + right modality (all-atom, metal) +
tractability, **not a singular SOTA claim**. Note every diffusion backbone hands
sequence design to **(Ligand)MPNN** — a smaller, more RL-tractable model that is the
cheapest first post-training target (and where the BoltzGen→LigandMPNN combo work
already showed coordination quality improving).

### Two post-training axes

1. **Sharpen the reward** — fine-tune the MLIP (MACE / UMA) on CSD-derived,
   DFT-labelled metal–organic geometries → a metal-coordination-accurate reward
   model (directly fixes the labile-donor weakness both MLIPs showed). CSD geometry →
   DFT single-point → fine-tune set; the standard universal-MLIP fine-tuning recipe.
2. **RL the generator** — use the verifier reward to fine-tune BoltzGen / LigandMPNN.
   RL fine-tuning of discrete protein diffusion + sequence models is established
   ([reward optimization with KL-to-pretrained for naturalness](https://arxiv.org/html/2410.13643v1);
   DDPP / RTB / DRAKES; DPO; [best-of-K diffusion alignment](https://arxiv.org/pdf/2501.15631)).
   Tractable order: **best-of-N / rejection sampling → DPO from verifier-ranked pairs
   → full policy-gradient RL** (full diffusion-backbone RL is the heaviest, do it last).

### Why the stack fits RLVR

Reward-hacking is the failure mode of single-reward RL. The
**2-independent-methods × 4-stages** design is the mitigation — a design must satisfy
geometry **and** bond-valence **and** physics **and** co-fold at once, far harder to
game than one z-score. The defense-in-depth that makes the verifier trustworthy is
exactly what makes it a robust RL reward.

## Structural-verifier vertical

MD → homogenise → FEM for mechanical soundness of periodic protein materials, sharing
the MLIP/MD substrate (see PR #5's scoping doc; gated on the RVE/periodicity
requirement).

## Where it lives

| PR | status | scope |
| --- | --- | --- |
| [#1](https://github.com/charleneleong-ai/ai4science/pull/1) | merged | core verifier stack — geometry oracle, generators, pipeline |
| [#2](https://github.com/charleneleong-ai/ai4science/pull/2) | open | `BoltzGenAdapter` — 2nd generator + mmCIF support |
| [#3](https://github.com/charleneleong-ai/ai4science/pull/3) | open | multi-metal selectivity (Ni/Cu/Co) |
| [#4](https://github.com/charleneleong-ai/ai4science/pull/4) | merged | reference-data fix (gitignore / package data) |
| [#5](https://github.com/charleneleong-ai/ai4science/pull/5) | open | MLIP physics tier + bond-valence + MLIP-MD |
| [#6](https://github.com/charleneleong-ai/ai4science/pull/6) | open | co-fold cross-check — Chai-1 + AllMetal3D |
| [#7](https://github.com/charleneleong-ai/ai4science/pull/7) | open | `CSDReference` (CSD/Mogul drop-in) |
| [#8](https://github.com/charleneleong-ai/ai4science/pull/8) | open | this spec — verifier-stack architecture + RLVR |

Stacked chain **#2 → #3 → #5 → #6** (merge bottom-up); **#7** and **#8** are
independent off `main`.
