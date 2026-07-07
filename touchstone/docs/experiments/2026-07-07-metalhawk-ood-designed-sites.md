# MetalHawk is confidently-OOD on designed sites — demote the tier (Ni, 2026-07-07)

End-to-end evaluation of the MetalHawk geometry-distortion tier (PR #51) on real BoltzGen
designs. Verdict: **demote it.** MetalHawk is a *learned* classifier and it extrapolates
confidently off its training manifold — the opposite of what a verifier needs. The analytic
[`coord_geometry`](../../src/touchstone/geometry/coordination.py) tier (polyhedron-RMSD) is
the geometry oracle for designed sites; MetalHawk stays as a licence-free, opt-in cross-check
for in-distribution (natural-like) sites only.

## What MetalHawk is meant to be
The open, learned stand-in for **Mogul** (the CSD-licensed geometry validator). Mogul is an
*empirical lookup*: it matches each geometric parameter against the distribution of that
fragment across the CSD and — critically — reports "insufficient data" when it has no match.
MetalHawk instead maps coordinates → one discrete geometry class (LIN/TRI/TET/SPL/SQP/TBP/OCT)
+ an entropy. That design choice is the whole story below: **a lookup abstains off-manifold;
a classifier hallucinates a confident class.**

## Setup
- **Env:** MetalHawk (github.com/vrettasm/MetalHawk) + pymol-open-source + its scikit-learn
  models, on `pi-a100-80gb`. Untangled a 3-layer ABI break — conda ships pymol compiled
  against numpy ≥1.23 (API 0xf), but MetalHawk pins numpy 1.21.5 (0xe), segfaulting pymol's
  `get_coords`; fixed with `numpy==1.23.5` + `scipy==1.9.3`.
- **Pool:** the round-4 fine-tuned Ni pool, 96 designs (`ni_motif`), gemmi-normalized
  CIF→PDB (BoltzGen CIFs crash pymol's metal parsing otherwise).
- **Pipeline:** `extract_metal_sites.py` (10 Å sphere) → `metalhawk.py` (classify) → parse,
  then [`verify_structure`](../../src/touchstone/service.py) with the predictions wired via
  [`metalhawk_score.py`](../../scripts/metalhawk_score.py) → `load_predictions` → `score_provider`.

## Result — confidently wrong
| model | class distribution (96 designs) | mean confidence |
|---|---|---|
| MetalPDB (`--no-csd`) | 95 OCT (CN6) · 1 TET (CN4) | 0.99 |
| CSD (`--csd`) | **96 OCT (CN6)** | ~1.00 (mean entropy 0.009) |

The designs are **physically CN4**. Independent distance check on samples: clean His/His/Asp/Cys
sites, all four donors ≤2.6 Å, **nothing 5th/6th within 3.0 Å**. MetalHawk isn't catching hidden
coordination — both its models call these tetrahedral sites octahedral at ~100% confidence.
It also warned on **all 96** that the sphere was "anomalously dense," consistent with the inputs
sitting off its CSD/MetalPDB manifold. CN agreement with the structure's own coordination: **1/96.**

## Consensus impact — it collapses the verdict
Two problems, in sequence. First the tier as written (`weak` on any CN mismatch) injected noise:
it downgraded 43 geometry-TRUST designs to `weak` on a spurious CN6-vs-CN4 mismatch. Adding an
**OOD gate** ([`MetalHawkVerifier`](../../src/touchstone/geometry/metalhawk.py): confident yet
Δcn ≥ 2 vs the physical shell ⇒ `defer`, not `weak`) made it *per-tier honest* — 95 `weak` → 95
`defer`. But `defer` collapses the consensus (defense-in-depth), so wiring the honest tier is
*worse* for the overall verdict:

| consensus | without MetalHawk | with MetalHawk wired |
|---|---|---|
| trust | 36 | **1** |
| weak | 51 | 0 |
| defer | 9 | **95** |

Whether it emits `weak` (noise) or `defer` (collapse), MetalHawk is **net-harmful as a consensus
voter on de-novo designs** — because it is OOD on ~every one of them.

## Decision — demote
- **`coord_geometry` is the geometry oracle for designs.** Its own docstring calls it "a gRMSD
  proxy, à la CheckMyMetal's geometry parameter" — it does MetalHawk's job analytically, with a
  real OOD gate (`ood_deg`), and degrades gracefully instead of hallucinating.
- **MetalHawk stays opt-in + experimental**, marked as such in the tier docstring and the service
  `_NEEDS_INPUT` description. Only wire it on in-distribution (natural-like) sites.
- **Kept from this work** (the tier is honest + runnable if used): the OOD gate; the fixed
  [`metalhawk_score.py`](../../scripts/metalhawk_score.py) (real `--inputdir/--outputdir` CLI +
  `file/prediction/entropy` parse + geometry→CN map; the original assumed a CLI/schema that don't
  exist); and `load_predictions`, which closes the JSON→verifier loop.

## Caveats
- **No Mogul head-to-head** — no CSD licence here; the comparison is by method (lookup-abstains
  vs classifier-extrapolates), not a benchmark.
- **One design class** (de-novo `ni_motif` Ni). MetalHawk may well be fine on natural sites; the
  claim is scoped to designed/OOD inputs, which is exactly touchstone's domain.
