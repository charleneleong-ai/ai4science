# BoltzGen for metal-binder design: scaffold a motif, don't "bind an ion"

**2026-06-25** · BoltzGen's touchstone-trusted Ni²⁺-site rate went from **0/16 → 3–5/12** by
changing only the *design spec* — giving it a coordination motif to scaffold around instead
of a bare ion. Nothing about the model or the verifier changed.

## The question

A first, naive comparison made BoltzGen look weak at metal-binder design:

| generator | spec | touchstone-trusted Ni sites |
| --- | --- | --- |
| RFdiffusionAA → LigandMPNN | de novo | 5 / 96 |
| BoltzGen | bare Ni²⁺ ion (`ligand: {ccd: NI}`) | **0 / 16** |

Every BoltzGen design came back `weak`/`defer` — CN 1–3, loose coordination. But BoltzGen is the
sponsor's tool, and that result is too convenient. *Why* did it fail?

## Diagnosis

The comparison wasn't apples-to-apples — it was posed to suit LigandMPNN and starve BoltzGen.

- **LigandMPNN** is inverse folding with *explicit* metal context: it picks residues + rotamers
  whose job is to coordinate the metal. Tight coordination is its training objective, so of course
  its output sits at CN 4–6.
- **BoltzGen** designs a binder *to a target* by reasoning about that target's **surface/interface**.
  A bare Ni²⁺ ion is a single, essentially surfaceless point — a near-degenerate target with nothing
  to grip. We handed it the one task LigandMPNN is purpose-built for, in the format worst-suited to it.

So "0/16" is a statement about the *spec*, not about BoltzGen.

## The fix — a theozyme spec

Pose metal design as BoltzGen's actual strength: scaffolding a fold around a specified
coordination motif. The self-contained way (no target structure file) is to fix the coordinating
residues in the **sequence string** and bond their donor atoms to the metal —
`boltzgen_in/ni_motif.yaml` (lives on the GPU box):

```yaml
entities:
  - protein:
      id: A
      sequence: "11H14H13D13C5"   # His12 His27 Asp41 Cys55, the rest designed
  - ligand:
      id: B
      ccd: NI
constraints:
  - bond: { atom1: [A, 12, NE2], atom2: [B, 1, NI] }   # His Nε2 → Ni
  - bond: { atom1: [A, 27, NE2], atom2: [B, 1, NI] }   # His Nε2 → Ni
  - bond: { atom1: [A, 41, OD1], atom2: [B, 1, NI] }   # Asp Oδ1 → Ni
  - bond: { atom1: [A, 55, SG],  atom2: [B, 1, NI] }   # Cys Sγ  → Ni
```

**Gotcha that cost an iteration:** fix the coordinating residues *in the sequence string*, not via
`residue_constraints`. The latter only constrain identity at design time — the sidechain atoms
aren't materialised, so a `bond` to `NE2`/`OD1`/`SG` fails at parse with `KeyError: ('A', 11, 'NE2')`.
Fixing them in the sequence makes the atoms exist. Validate with `boltzgen check <yaml>` before running.

## Result

| BoltzGen run | spec | trusted | donors on trusted designs |
| --- | --- | --- | --- |
| bare ion | wrong | 0 / 16 | — |
| **+ theozyme** | right | **5 / 12** design step, **3 / 12** refold | `NNOS`, `NNOOS` (His·His·Asp·Cys), CN 4–5 |

The trusted designs coordinate the Ni with exactly the His/His/Asp/Cys donors specified, at the
geometry of real Ni²⁺ sites. The refold step is slightly more conservative than the raw design step
(3 vs 5) — re-prediction loosens some sites, which touchstone honestly catches. Three runs are logged
to W&B for side-by-side comparison (`ligmpnn-nickel-filter`, `boltzgen-nickel-filter`, `boltzgen-motif`).

## Enhancing it: hand the scaffold to LigandMPNN — but keep the motif

The obvious refinement is to pass BoltzGen's designs through LigandMPNN (which is metal-aware).
It only helps if done in *refine* mode, not *redesign* mode:

| approach | trusted |
| --- | --- |
| BoltzGen theozyme alone | 3 / 12 (25%) |
| → LigandMPNN, full redesign | 2 / 24 (8%) — **discards the motif** |
| → LigandMPNN, coordinators fixed | **8 / 24 (33%, 0 defer)** |

Unlike RFdiffusionAA — which emits bare backbones, so LigandMPNN *adds* the coordinating
sidechains (DEFER → trust) — BoltzGen already places them. A full LigandMPNN redesign therefore
*throws away* BoltzGen's His/His/Asp/Cys and coordinates worse. Fixing the four coordinating
residues (`--fixed_residues "A12 A27 A41 A55"`) and letting LigandMPNN only re-pack their rotamers
+ optimise the scaffold tightens the geometry and removes the failures (0 defer). The lesson:
pair a *full-design* generator with a packer in **refine** mode, not **redesign** mode — and the
verifier is what tells the two apart.

## Why it matters

The verifier did more than score — it **diagnosed the misuse and confirmed the fix**. That is the
verifier-first thesis doing real work: an independent verifier doesn't just gatekeep a generator, it
tells you *how to use the generator correctly*, then proves the change worked. The same
[`GeometryVerifier`](../src/touchstone/geometry/verifier.py) judged all three runs through the
identical [`BoltzGenAdapter`](../src/touchstone/generators.py) path — no special-casing.

For the sponsor's framing: **BoltzGen designs trustworthy metalloproteins when posed as
motif-scaffolding**, and touchstone is the independent check that proves it.

## Reproduce

```bash
# on the GPU box (conda env bg)
boltzgen check boltzgen_in/ni_motif.yaml
boltzgen run boltzgen_in/ni_motif.yaml --output boltzgen_motif_out \
  --protocol protein-anything --num_designs 12 --steps design inverse_folding folding

# verify + log (touchstone, format-agnostic)
uv run --extra viz python scripts/log_wandb.py <refold_cif_dir> --generator boltzgen-motif
```
