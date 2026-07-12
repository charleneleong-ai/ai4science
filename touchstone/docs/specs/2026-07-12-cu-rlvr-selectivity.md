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

## Open items before a real run

- **A Cu-native motif spec** for generation (N-rich square-planar / type-1), not the `ni_motif`
  theozyme. The reward can only select from what the generator proposes.
- **Audit the Cu²⁺ geometry prior.** The active CSD reference reports Cu²⁺ modal **CN = 3** (vs Ni/Co
  at 6), which looks wrong — Cu²⁺ is typically 4–6. A bad prior blunts the geometry tier for Cu;
  rebuild from MetalPDB ([`scripts/build_metalpdb_reference.py`](../../scripts/build_metalpdb_reference.py)).
- **Jahn-Teller-aware ideal geometry** — [`coord_geometry`](../../src/touchstone/geometry/coordination.py)
  scores against symmetric ideal polyhedra; a JT-elongated octahedron would judge Cu fairly.
