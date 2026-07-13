# Cu RLVR — a selectivity-aware reward (2026-07-12)

> ## ⚠️ SUPERSEDED — the reward this spec is built on does not work
>
> The MLIP metal-swap selectivity tier **cannot rank metals**: MACE-MP and OrbMol both fail the
> Irving–Williams series (neither puts Cu²⁺ at the peak), even given correct charge, per-metal spin
> and full relaxation. Metal preference is a ligand-field effect these potentials do not carry.
> The tier is now **gated off** and `defer`s, so `--selectivity` is inert — see
> [**MLIPs cannot rank the divalent 3d series**](../experiments/2026-07-13-mlip-cannot-rank-metals.md).
>
> **The plan below is not executable as written.** The reward needs a backbone with ligand-field
> physics (DFT with explicit spin states), not an MLIP. The Cu-native motif reasoning still stands;
> the *selection signal* does not.

Scaffold for applying the proven Ni RLVR loop
([`docs/experiments/2026-07-01-rlvr-boltzgen-round1.md`](../experiments/2026-07-01-rlvr-boltzgen-round1.md))
to **Cu²⁺** — with the one lever the Ni rounds lacked: **selectivity in the reward**.

## Why Cu isn't just "Ni with a different label"

Verifying a Ni-designed site (`ni_motif_37`) as Cu²⁺ under the full stack exposed the problem:

| tier | Ni²⁺ (as designed) | Cu²⁺ (same coords) |
|---|---|---|
| bond_valence | 1.83 | 2.01 |
| precedent | 4× Ni-N2O1S1 | 3× Cu-N2O1S1 |
| **mlip** | site held, ΔE_bind **−4.19 eV** | lost a donor, ΔE_bind **−4.86 eV** |
| mlip_md | 94% | 65% |
| **selectivity** | **ΔE favours Co²⁺ (margin −0.64 eV)** — the Ni target is out-competed | — |

Cu binds this site **more strongly than the Ni it was designed for** (−4.86 vs −4.19 eV) yet is less
dynamically stable. That's the whole problem in one line: **the geometry/physics tiers reward a site
that binds — they don't reward a site that binds the *target* metal preferentially.** In a mixed
leachate (Ni/Cu/Co/Fe), a strong-but-unselective binder is useless.

Cu²⁺ also has genuinely different coordination chemistry: it favours N-rich square-planar (His₃/His₄)
or type-1 "blue copper" (His₂Cys/Met) over the octahedral His/His/Asp/Cys `ni_motif`, and its
octahedra are Jahn-Teller-elongated. So Cu needs a **Cu-native motif spec**, not a relabelled Ni one.

## The reward

Same shape as the Ni reward (mean verifier score × consensus weight, see
[`reward_from_result`](../../src/touchstone/reward.py)) — but with the **metal-swap ΔE selectivity
tier** ([`MLIPSelectivityVerifier`](../../src/touchstone/physics/selectivity.py)) folded in, so the
consensus/reward only rewards designs whose **target metal binds most favourably** among competitors:

```bash
# in the mace env on the GPU box
uv run python scripts/rlvr_select.py \
    --npz-dir <fold_out_npz> --cif-dir <refold_cif> --out cu_round1 \
    --metal Cu2+ --deep --selectivity Ni2+,Cu2+,Co2+ --keep trust
```

[`rank_structures`](../../src/touchstone/reward.py) now threads `selectivity_metals` through to
`verify_structure`, so `selectivity` joins the consensus: a design where Ni or Co out-competes Cu
scores `weak` on that tier and is demoted — exactly the discrimination geometry can't make.
`--selectivity` requires `--deep` (it's an MLIP tier).

## The loop

Unchanged from the Ni arc — only the reward's metal + the selectivity tier differ:

    BoltzGen generate (Cu-native motif)  →  rlvr_select --metal Cu2+ --selectivity …  →
    winners_to_targets  →  BoltzGen resume-train  →  re-verify a fresh pool  →  repeat

## A Cu-native motif is required

The reward can only select from what the generator proposes, and the `ni_motif` theozyme
(octahedral His/His/Asp/Cys) is a *Ni* site. A Cu spec should target Cu²⁺'s actual preferences:

- **Type-2 / square-planar**: N-rich (His₃, His₄, or His₂ + Asp/Glu) — the common Cu²⁺ protein site.
- **Type-1 "blue copper"**: His₂Cys + axial Met/Gln — high covalency, strongly Cu-preferring, and
  the donor set most likely to *discriminate* Cu from Ni/Co (soft Cys/Met thiolate/thioether).

The type-1 donor set is the one to try first: selectivity is a **donor-identity (HSAB)** effect, so a
soft-donor motif is what gives the metal-swap ΔE tier something to reward.

## Audit of the Cu²⁺ geometry prior (done — two real findings)

| prior | Cu²⁺ modal CN | bond mean ± std | Ni²⁺ std (for scale) |
|---|---|---|---|
| **CSD** (active) | **3** | 2.118 ± **0.211** | 0.108 |
| PDB | 4 | 2.144 ± 0.179 | 0.183 |

1. ~~**Cu's bond-length std is ~2× Ni's — that is the Jahn-Teller signature.**~~ **Retracted — this was
   an artifact of the CSD prior, and it dissolves once the domain is fixed.**

   JT elongation is a **six-coordinate** effect (≈4 short + 2 long bonds). In the *protein* domain
   only **3 of 309 Cu²⁺ sites (1.0%)** are 6-coordinate — protein copper is 3–4 coordinate
   (type-1/type-2), and 99% of it structurally *cannot* Jahn-Teller distort. The physics is real
   where it can occur (in those 3 sites: 4 short bonds at 2.107 Å, 2 long at 2.630 Å — a textbook
   **+0.52 Å** elongation), but it is not what makes the Cu prior wide.

   The width (±0.216 Å) is **donor heterogeneity**: Cu–N 2.119, Cu–O 2.230, Cu–S 2.273 Å. The
   distribution is broad-unimodal, not bimodal. CSD is full of octahedral Cu chelates, so the JT
   reading was correct *there* — and irrelevant here. **Modelling Cu bimodally would be fixing a
   small-molecule artifact inside a protein prior.** Item closed, not deferred.
2. **The active prior is the wrong domain.** [`best_reference`](../../src/touchstone/geometry/reference.py)
   already prefers MetalPDB → CSD → PDB, but `metalpdb_reference.json` was never built/bundled — so
   touchstone verifies **protein** metal sites against a **small-molecule crystal** (CSD) prior. That
   also explains the odd Cu modal CN=3: CSD Cu²⁺ is full of low-coordinate chelates, whereas protein
   Cu²⁺ is typically 4–5.

   ✅ **Fixed and worse than suspected.** The MetalPDB prior is built and is now the default. Scored
   against real metalloprotein sites, the old CSD prior trusted only **41% of real Ni²⁺ sites** and
   deferred **19.5%** of them as off-manifold; the corrected prior trusts 96.6%. It was centred
   0.116 Å short — over 1σ of its own std — so the RLVR reward had been steering BoltzGen toward
   small-molecule bond lengths. See
   [**the geometry prior was measuring the wrong domain**](../experiments/2026-07-13-geometry-prior-wrong-domain.md).

   **This bites Cu hardest, which is not obvious from the trust rates.** CSD trusts 97% of real Cu²⁺
   sites, so it *looks* fine — but only because it is too blunt to reject anything: its `cn_range`
   runs to 6, while protein Cu²⁺ is 3–5 coordinate. Given an **octahedral CN6 copper** site (aqua
   geometry; no protein Cu site is octahedral) the CSD prior returns `plausible (0.0σ)` — a perfect
   score. So **the Cu geometry reward on the old prior would have rewarded octahedral copper**, when
   the whole point of the Cu motif work below is square-planar / type-1 geometry. The corrected prior
   rejects it (`coordination number outside observed range`).

**Both items are now closed** — one fixed, one retracted. The prior audit is done.
