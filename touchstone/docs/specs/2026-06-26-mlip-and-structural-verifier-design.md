# MLIP tier + structural-verifier vertical — design

**Date:** 2026-06-26
**Status:** scoping (phase 1 in progress)
**Builds on:** [`2026-06-24-touchstone-design.md`](2026-06-24-touchstone-design.md)

## Why

Touchstone's verdict today rests on a geometry oracle, an independent co-fold
(Boltz-2), and a semi-empirical physics tie-breaker (`xtb` GFN2). The physics
tier is the weakest link: GFN-FF segfaulted repeatedly, GFN2 is a coarse
approximation of a transition-metal coordination shell, and neither spans the
materials regime. The goal is to **enhance the verification layer so the verdict
correlates better with wet-lab outcome** — sharper physics, and a path to a
second (mechanical) verification axis.

Honest framing of "verify wet-lab success": a physics tier sharpens the
**thermodynamic binding + geometric** signal only. Expression, fold stability,
solubility, kinetics, and the assay are invisible to it. So the deliverable is a
**better-calibrated triage signal** (spend scarce wet-lab slots on the right
designs), not a success oracle. Encouragingly, physics/confidence metrics *do*
correlate with wet-lab hits across large designed-binder sets — but the
correlation is only *proven* once the tier is calibrated against real outcomes,
which touchstone has none of yet. Calibration is the eventual unlock, tracked as
phase 3.

## Architecture: two verticals, one substrate

An MLIP (machine-learned interatomic potential — DFT-accuracy energies/forces at
MD speed) is the shared engine that unlocks both a sharper binding check and an
entirely new mechanical check:

```
                           ┌─ binding verifier:   geometry → cofold → MLIP energy    "does it grab the metal?"
design ─→ MLIP-MD ─→ ┤
                           └─ structural verifier: MD → homogenise (Cᵢⱼ) → FEM        "does the material hold load?"
```

Both stay generator-blind. The MLIP tier is not just an upgrade to the binding
check — it is the substrate that makes the structural axis possible at all.

## Phase 1 — MLIP binding tier (in progress)

Drop-in alongside the `xtb` step, reusing the existing cluster extraction
([`scripts/extract_cluster.py`](../../scripts/extract_cluster.py)): relax the
metal + first/second shell under an MLIP and score the metal site.

**What it scores (both):**
- **Stability gate** — relax the coordination cluster; flag designs whose metal
  site drifts or loses a donor. Drift is measured over the *metal + first-shell*
  atoms only; cluster-edge capping atoms always wander and are excluded.
- **Interaction-energy ranking** — ΔE of binding the metal (complex vs apo +
  ion). Ranking-only: absolute solution-phase binding energies from MLIPs are
  unreliable, so this orders designs rather than gating them, pending
  calibration.

**Model choice — benchmark before wiring.** Candidates: MACE-MP (Materials
Project foundation, permissive, has transition metals), MACE-OFF (organic-only —
no Ni, so it cannot do a metalloprotein site alone), UMA (FAIR universal, spans
organic + metal, gated download). Benchmark = relax the Ni cluster and measure
agreement with the `xtb` GFN2 optimum (`xtb_work/cluster_opt.pdb`).

Findings so far (2026-06-26, A100, `mlip` conda env):
- **MACE-MP** reproduces the GFN2 metal-site geometry to **0.056 Å MAE on the
  three tight Ni–donor bonds in ~8 s**, but **disagrees with xtb on the labile
  4th donor**: GFN2 *retains* it, tightening it to 2.32 Å (CN 4); MACE-MP
  *expels* it past 2.8 Å (CN 4→3). A genuine accuracy gap on labile metal–organic
  coordination — the bond type CN-based verdicts hinge on — not a harness
  artifact. (Which model is right is unknown without DFT; GFN2 is itself
  approximate, so the disagreement is arguably useful "borderline coordination"
  signal.) Permissive, zero gating.
- **UMA vs MACE-MP cannot share a conda env** — fairchem pulls `e3nn 0.6` /
  `torch 2.8`, mace pins `e3nn 0.4.4`. Two envs required; a deployment cost that
  favours MACE-MP unless UMA is clearly more accurate on the labile donor.
- **UMA** (gated; `facebook/UMA`, FAIR Chemistry License w/ acceptable-use +
  geographic restrictions — separate grant from the CC-BY-4.0 OMAT24 dataset)
  reproduces the tight donors to **0.083 Å in 16 s**, and **independently makes
  the same CN 4→3 call** — it expels the labile 4th donor too. Two foundation
  MLIPs agree it is not bound; only semi-empirical xtb GFN2 retains it, so xtb is
  the likely outlier (over-binding), and CN 4→3 is a defensible "borderline
  coordination" signal rather than a MACE artifact.

**Decision: MACE-MP-0** as the backbone — more accurate on the tight donors
(0.056 vs 0.083 Å), 2× faster, and MIT vs UMA's gated/AUP license. UMA stays the
swappable independent check (the pluggable backbone), which is how it earned its
keep here: corroborating the labile-donor verdict rather than overturning it.

**Integration:** runs via ASE (`mace-torch` / `fairchem`) on GPU, **lazy-imported**
like `wandb` so the numpy-only core stays light; gated behind a `touchstone[mlip]`
extra. Emits a sub-verdict of the same shape the stack already consumes.

## Phase 2 — structural-verifier vertical (spike, deferred)

`design → MLIP-MD → homogenise → FEM` — a validated hierarchical (FE²)
multiscale workflow (proven <5% vs experiment for lattice/TPMS bone scaffolds).
Each arrow:
- **MLIP-MD** — thermal + steered MD of the designed material.
- **Homogenise → Cᵢⱼ (E, ν)** — apply affine strains, read stress off the virial.
  `matscipy.elasticity` does this from an MLIP in a few lines.
- **FEM** — feed Cᵢⱼ into FEniCS/dolfinx; macro-scale scaffold load tests.

**Hard prerequisite — periodicity (RVE).** A single binder is not a continuum and
has no meaningful modulus. Homogenisation needs a *repeating* representative
volume, so this axis applies only to a **self-assembling / periodic protein
material** (designed lattice, fiber, hydrogel). BoltzGen designs binders, not
structural lattices — so the front-end generator is wrong as drawn. Resolution:
a protein-lattice/oligomer/fiber generator builds the scaffold; a binder model
(BoltzGen) *decorates* it with the metal-binding function.

**Other caveats:**
- **Hydration dominates** — protein moduli differ by orders of magnitude dry vs
  solvated; must model physiological hydration.
- **Validation gap** — predicted modulus needs AFM / nanoindentation / rheology
  to calibrate; protein-specific MLIP elastic-constant accuracy is unproven
  (benchmarks are inorganic crystals).
- **Cost** — converged Cᵢⱼ needs ns-scale MD across multiple strain states.

**Demo vs program:** a contained spike — one *periodic* designed protein material
→ MLIP-MD → Cᵢⱼ → a single FEM load test — is a striking, doable proof-of-concept.
The robust general pipeline is a research program, not a weekend build.

## Sequencing

1. **Phase 1** — finish the MLIP binding tier (shared substrate; immediate
   triage win). ← current
2. **Phase 2** — structural-verifier spike on a periodic material.
3. **Phase 3** — calibrate either axis against wet-lab / experimental data.

## Out of scope (for now)

- Continuum corrosion / degradation modelling (phase-field) — relevant to the
  separate biodegradable-metals framing, not to binding or mechanical soundness.
- Electronic conductivity (quantum transport) — needed for the bioelectronics
  framing; neither geometry nor an MLIP captures it.
