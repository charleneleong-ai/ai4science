# touchstone — design spec

**Date:** 2026-06-24 · **Status:** approved-pending-review · **Author:** Charlene

> *A touchstone was historically used to assay precious metals — the independent reference
> standard against which a sample is judged. That is exactly this component's job.*

## Summary

A **generator-agnostic verifier** for designed metal-binding proteins: any generator
produces a candidate binder; an independent `GeometryVerifier` scores its metal-coordination
geometry against a reference distribution, flags out-of-distribution (extreme-leachate)
inputs, and returns a ranked list with a **trust / defer** verdict per design.

Domain: protein–metal binders for **critical-mineral recovery**. A standalone pet project
that makes the verifier-first thesis (see [`../../README.md`](../../README.md)) concrete.

## Motivation

When the generator is commoditised — every effort can reach for the same de-novo binder
model (RFdiffusionAA, BoltzGen) — the only axis left to differentiate on is the **verifier**.
The easy move is to read the generator's own confidence; the contrarian, defensible move is an
*independent* geometry oracle that knows when a designed site is physically implausible or
when the input has gone out-of-distribution (hot, acidic, saline leachate). That is the whole
thesis of this repo, reduced to one shippable component.

## Scope

**In scope** — exactly this pipeline, nothing more:

```
Generator (adapter)                    Verifier (the asset, generator-blind)
  BoltzGenAdapter   ──┐                  GeometryVerifier
  RFdiffusionAdapter ─┼─▶ BinderDesign ─▶   ↳ score metal-site geometry vs reference dist
  MockGenerator     ──┘   (seq, struct,     ↳ OOD flag: distance under extreme-condition shift
                           target_metal)  ──▶ ranked, with "trust / defer" verdict
```

**Out of scope (YAGNI):** no multi-domain framework (no `conformal`/`rlvr` modules — those
stay as README context only); no wet-lab step; no generator *training*; no novel metal
chemistry. The broader 5-project pattern is positioning, not code.

## Architecture

Three small units with one contract between the two halves.

### `BinderDesign` (the contract)
Dataclass the generator emits and the verifier consumes — the only coupling point.

| field | meaning |
| --- | --- |
| `sequence` | designed protein sequence |
| `structure` | 3D coords (PDB / mmCIF), incl. the bound metal |
| `target_metal` | element + oxidation state (e.g. `Cu2+`) |
| `generator_confidence` | the generator's *own* score — recorded, never trusted by the verifier |

### `Generator` (swappable, generator-agnostic)
Protocol: `design(target) -> list[BinderDesign]`. Three implementations:
- `RFdiffusionAdapter` — **primary POC generator**, runs on the A100. Purpose-built for
  metal-coordination design: it places the metal atom explicitly (exactly what the verifier
  needs to score a site) and ships a bundled nickel example for an instant first run.
- `BoltzGenAdapter` — a second real generator. Swapping it in *proves* the verifier is
  generator-blind (nothing downstream changes).
- `MockGenerator` — instant, deterministic; lets the verifier be built/tested without GPU.

### `GeometryVerifier` (the asset)
`verify(design) -> Verdict`. Independent of how the design was made. Steps:
1. **Extract the coordination site** — metal centre, coordinating atoms, bond lengths, angles,
   coordination number.
2. **Score vs reference distribution** — how typical is this geometry for `target_metal`?
   Reference is *pluggable*: a public PDB pull now, swappable for **CSD/Mogul**
   (real empirical bond-length/angle distributions). Same interface either way.
3. **OOD flag** — perturb the input toward extreme leachate (pH / temperature / salinity proxy)
   and measure how far the site moves off the reference manifold. Large shift ⇒ defer.

`Verdict` = `score: float`, `trust: bool`, `ood: bool`, `reason: str`. The pipeline returns
designs ranked by `score`, each with its verdict.

## Reference data (pluggable oracle)

- **PDB pull:** scripted pull of **Ni²⁺/Cu²⁺** coordination sites from the public PDB →
  empirical bond-length/coordination-number distributions. Real logic, no license needed.
- **CSD/Mogul:** swap the reference provider for the **CSD Python API + Mogul**. The
  `GeometryVerifier` does not change — only the `ReferenceDistribution` implementation behind it.

**Metal:** `Ni2+` for the **POC** — RFdiffusion ships a nickel design example, so the generator
runs end-to-end with zero config. `Cu2+` stays an option for the verifier's reference data
(richer PDB + CSD coverage); both are real critical-mineral recovery targets.

## Environment

- Host: `pi-a100-80gb` (216.81.200.11) — A100 80GB, idle, miniconda present. **Disk at 90%
  (~80 GB free)** — watch it when pulling BoltzGen weights.
- Python 3.11+, `uv` for the `touchstone` package; conda env for the GPU/BoltzGen side.

## Setup sequence (run-for-real)

1. **A100 env + run RFdiffusionAA's bundled nickel example end-to-end** — confirms GPU + I/O,
   produces a real Ni²⁺-coordinating binder with the metal atom placed (drives what
   `RFdiffusionAdapter` emits). BoltzGen env is a later step.
2. **`Generator` interface + `RFdiffusionAdapter` + `MockGenerator`.**
3. **`GeometryVerifier` + the PDB-Cu²⁺ `ReferenceDistribution` stand-in + OOD flag.**
4. **One toy target end-to-end** — generator → verifier → ranked list with trust verdicts.
   Then: swap mock→real CSD reference and toy→real target.

## Testing

- `MockGenerator` + a hand-built `BinderDesign` with a *known-good* and a *known-bad* Cu²⁺
  geometry → verifier ranks good above bad, flags the distorted one.
- OOD: an in-distribution site scores `ood=False`; the same site under an extreme-condition
  perturbation scores `ood=True`. Assert the *crossover*, not just a smoke run.
- Reference provider swap (PDB stand-in ↔ a second mock provider) leaves the verifier API
  unchanged — proves the oracle is pluggable.
- Generator swap (`MockGenerator` ↔ `BoltzGenAdapter` interface) leaves the verifier unchanged
  — proves generator-blindness.

## Success criteria

- The generator (RFdiffusionAA) runs on the A100 and its I/O is understood.
- `generator → BinderDesign → GeometryVerifier → ranked trust/defer` runs end-to-end on the
  PDB pull, generator-agnostic, with the tests green.
- Swapping in the real CSD/Mogul reference and a real target is a config change, not a rewrite.
