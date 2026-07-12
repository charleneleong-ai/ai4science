# Cu RLVR — a selectivity-aware reward (2026-07-12)

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

1. **Cu's bond-length std is ~2× Ni's — that is the Jahn-Teller signature.** A JT-elongated Cu²⁺
   octahedron has *bimodal* bonds (≈4 short + 2 long), and the prior crushes them into one Gaussian.
   The result isn't a wrong mean so much as a **blunt** one: the std is so wide that almost any Cu
   bond length passes the z-score, so the geometry tier barely discriminates for Cu.
   *Correction to an earlier note in this doc:* the fix is **not** a JT-elongated ideal in
   [`coord_geometry`](../../src/touchstone/geometry/coordination.py) — JT changes bond **lengths**, not
   **angles** (a JT octahedron still has ~90°/180°), so the angular tier is already fine. The gap is
   the single-Gaussian **bond-length** model in the geometry tier.
2. **The active prior is the wrong domain.** [`best_reference`](../../src/touchstone/geometry/reference.py)
   already prefers MetalPDB → CSD → PDB, but `metalpdb_reference.json` was never built/bundled — so
   touchstone verifies **protein** metal sites against a **small-molecule crystal** (CSD) prior. That
   also explains the odd Cu modal CN=3: CSD Cu²⁺ is full of low-coordinate chelates, whereas protein
   Cu²⁺ is typically 4–5.

**Next step (its own PR — it changes the active prior for every metal and every verdict):** run
[`scripts/build_metalpdb_reference.py`](../../scripts/build_metalpdb_reference.py) and bundle
`metalpdb_reference.json`, giving Ni/Cu/Co metalloprotein-specific priors. Model the Cu bond length
bimodally (or widen only the axial component) rather than with one Gaussian.
