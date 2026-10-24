# Calibration: what the scores are worth, measured before using them

Two gates decide this entry — which designs get ranked, and whether any of them activates the
receptor. Neither can be trusted to a number that has not been checked, so both were calibrated
against known answers before any design was generated. Reproduce with
[`calibrate_scores.py`](../scripts/calibrate_scores.py) and
[`activation_geometry.py`](../scripts/activation_geometry.py).

## 1. Ranking: no available score is worth much, and ipTM is worse than nothing

The Proteinbase corpus publishes predicted scores *and* wet-lab outcomes for the same designs,
so the organisers' own conclusion that their metrics did not predict binding is checkable. It
checks out, and the detail matters. **4,069 designs, 2,214 with Boltz-2 scores, 3,169 with a
binding label, 27 targets.**

Pooled AUROC across targets is misleading — easy targets inflate it — so the column that counts
is **within-target**, weighted by the comparisons each target supports.

| score | binding, within-target | expression, within-target |
| --- | --- | --- |
| **boltz2_complex_iplddt** | **0.652** | 0.692 |
| boltz2_plddt | 0.628 | 0.710 |
| boltz2_complex_plddt | 0.607 | **0.729** |
| boltz2_iptm | 0.607 | 0.624 |
| boltz2_lis | 0.593 | 0.647 |
| esmfold_plddt | 0.578 | 0.672 |
| seqidentity | 0.537 | 0.321 |
| **glu_fraction** | 0.510 | **0.724** |
| lys_fraction | 0.510 | 0.605 |
| proteinmpnn_score | 0.480 | 0.322 |
| boltz2_complex_pde | 0.475 | 0.442 |

Three things follow.

**Nothing reaches 0.66 for binding.** The best available signal is interface pLDDT at 0.652.
Ranking 20 slots on that is defensible; ranking them on ipTM is not.

**Every score predicts expression better than binding.** `boltz2_complex_plddt` reaches 0.729
against expression but 0.607 against binding. So a good part of what these scores detect is
"this protein folds", not "this protein binds the target" — which is worth knowing, because
expression is a genuine gate that 9% of designs fail.

**The composition result replicates.** Glutamate fraction predicts expression at 0.724
within-target, against 0.77 published. (I first tested it against binding, where it is 0.510 —
that was my error, not a refutation. It is an expression predictor and a good one.) It costs
nothing to compute, so it goes in as a pre-filter.

### ipTM does not rank binders, and over most of its range it ranks them backwards

| ipTM band | n | observed binding rate |
| --- | --- | --- |
| 0.0–0.2 | 47 | **46.8%** |
| 0.2–0.4 | 151 | 40.4% |
| 0.4–0.6 | 171 | 34.5% |
| 0.6–0.7 | 169 | 34.9% |
| 0.7–0.8 | 207 | 33.8% |
| 0.8–0.9 | 482 | **21.2%** |
| 0.9–1.0 | 987 | 38.7% |

The rate *falls* from 46.8% to 21.2% as ipTM rises to 0.9, then jumps back. Against expression
the same bands are flat at 86–100%, so this is specific to binding.

**This is not a causal claim.** The corpus pools many campaigns and targets, band populations are
wildly uneven (987 designs at 0.9–1.0 against 47 at 0.0–0.2), and entry into the scored subset is
not random — designers plausibly submitted high-ipTM designs for the hardest targets, which would
produce exactly this shape. The defensible conclusion is narrower and sufficient: **a high ipTM
band does not contain a higher fraction of binders, so ranking by ipTM would not have helped.**

### A data-quality finding worth reporting upstream

`boltz2_pdockq` is **constant at 0.0183** across all 2,214 scored designs, and `boltz2_pdockq2`
constant at 0.0073. These are placeholders, not weak predictors — an AUROC of exactly 0.500 on
1.1M comparisons is the signature. `calibrate_scores.py` now detects and names constant fields
rather than reporting them as "no signal", because the two statements mean different things.
`seqidentity` is published as a nested object rather than a number and needs unwrapping.

## 2. Activation: the geometric gate orders real structures by function

No tool predicts agonism, so the gate is written explicitly — six pharmacophore features ranked
by the mutational cost of losing them (see [`strategy.md`](strategy.md) §2), plus insertion depth
below the vestibule rim. Cutoffs were set 1.0 Å looser than the native complex achieves, from
6X18 alone.

Scored on **experimental** structures of ligands whose function is known. Only 6X18 informed the
cutoffs; the rest were blind.

| structure | ligand | known function | depth Å | features | tier |
| --- | --- | --- | --- | --- | --- |
| 6X18 | GLP-1(7-36) | full agonist, EC50 0.9 nM | 20.6 | **6/6** | deep |
| 6B3J | exendin-P5 | full agonist, 0.85 nM, Emax 102% | 23.6 | 4/6 | deep |
| 7LLL | exendin-4 | full agonist, 0.09 nM | 18.9 | 4/6 | deep |
| 9IVG | GLP-1(9-36) | **weak partial, 363 nM** | 14.2 | **1/6** | shelf |

Plus a truncation ladder built from the native peptide in its own crystallographic pose, which
supplies the groove-only negative that no experimental structure provides (every published
GLP1R antagonist was solved on the isolated ECD, where there is no bundle to measure against):

| control | residues | depth Å | features |
| --- | --- | --- | --- |
| GLP-1 16–36 — groove only | 21 | 7.3 | **0/6** |
| GLP-1 13–36 | 24 | 11.7 | 1/6 |
| GLP-1 9–36, held in the agonist pose | 28 | 19.0 | 5/6 |

The three full agonists separate from the weak partial, and the groove-only fragment scores zero
features. The tier thresholds (deep ≥16 Å and ≥3 features; shelf ≥10 Å) are set from this ladder
and are not tuned afterwards.

The feature pattern also tracks published mechanism rather than just passing or failing.
Exendin-4 misses the Q234 hydrogen bond at 4.93 Å, matching its reported "lack of persistent
receptor interactions" for its first residue; exendin-P5 misses the R190 anion at 5.03 Å and is
the known G-protein-**biased** agonist. That the gate reproduces those distinctions without being
told about them is the strongest evidence it measures something real.

### The limitation, which is severe

**The gate scores the pose it is given, not whether the design would adopt it.** GLP-1(9-36) held
in the agonist pose scores 5/6 and "deep"; the same sequence in its *own* experimentally
determined pose scores 1/6 and "shelf", because it actually rebinds 7–12 Å further out. The gap
between those two rows is the entire phenomenon the gate is meant to detect, and the gate cannot
tell them apart from coordinates alone.

So a design that a predictor happens to place in an agonist-like pose will score well whether or
not it would really bind that way. Mitigations, in order of how much they are worth:

1. **Cross-model pose agreement** — the same design co-folded by independent predictors, scored
   on whether they agree about the pose, not on either one's confidence. `chai.sif` supplies
   Chai-1 as the independent model, so BoltzGen designs are never ranked only by the Boltz-2
   weights that made them.
2. **Many seeds.** Seed dependence is reported as substantial for every predictor on
   GPCR–peptide complexes, and the TNF-α campaign found three samples reversed four verdicts
   that twenty samples settled.
3. **Report the funnel with denominators**, including how many designs failed the gate, rather
   than only the survivors.

## 3. What this fixes in the plan

- Rank on **interface pLDDT with cross-model pose agreement**, never ipTM.
- Pre-filter on **glutamate fraction** for expression — replicated, free, and a real 9% gate.
- Spend nothing on Boltz-2's affinity head: it requires a small-molecule ligand chain, caps at
  128 atoms and was trained at ≤56, against ~450 for a 60-aa design.
- Treat the geometric gate as the only function-relevant signal available, and state plainly that
  it measures pose quality rather than pose probability.

## 4. Not established

- None of this is GLP1R-specific. The corpus has **zero GLP1R rows and no GPCR of any kind** —
  this is Adaptyv's first GPCR and first detergent-solubilised membrane target, so there is no
  target-matched calibration set and the AUROCs above are borrowed from soluble targets.
- The gate has never been tested on a *predicted* complex, only on experimental coordinates. Its
  behaviour on a Boltz-2 or Chai-1 model of a de novo design is unmeasured.
- No full-length experimental structure of a GLP1R antagonist exists, so the groove-only negative
  is a construct from the native pose rather than an observed antagonist complex.

## 5. The gate on *predicted* structures — the gap closed, with a bad answer

Section 2 validated the activation gate on experimental coordinates. The design pipeline will
only ever see predicted ones, so the same panel was co-folded with Boltz-2 against the 409-aa
assay construct, 5 diffusion samples each, and scored with the identical gate.

| ligand | known function | mean features | max features | deep / 5 |
| --- | --- | --- | --- | --- |
| exendin-P5 | full agonist | **3.8** | 4 | 5/5 |
| **GLP-1(15-36)** | **inactive** | **3.0** | 3 | **4/5** |
| exendin-4 | full agonist | 2.0 | 5 | 2/5 |
| GLP-1(9-36) | weak partial | 1.8 | 3 | 3/5 |
| GLP-1(7-36) | full agonist | 1.6 | 4 | 2/5 |
| exendin-4(9-39) | antagonist | 1.0 | 3 | 0/5 |
| GLP-1 H7A | inactive | **0.2** | 1 | 0/5 |

**The gate works as a negative filter and fails as a positive one.** The two clear negatives sit
at the bottom on mean features and score 0/5 deep. But an **inactive 15-mer ranks second**, above
two genuine full agonists, and the native agonist GLP-1(7-36) ranks fifth of seven.

The mechanism matters for design. GLP-1(15-36) is the ECD-binding half with the activation
segment removed; having no N-terminal half, it has nothing anchoring it to the groove, and
Boltz-2 slots it into the bundle instead. **Short designs are therefore at specific risk of
spurious deep placement** — the prediction, not the gate, is what fails.

This also reproduces the published warning independently: on ligands over 50 residues, 35.5% of
Boltz-2 predictions score DockQ < 0.23, and roughly a fifth of Boltz GPCR predictions land in the
wrong pocket with plausible confidence. Scoring `model_0` alone gave GLP-1(7-36) **0/6** against
its experimental **6/6**; only looking across all five samples recovered any signal at all,
which is the same single-sample error corrected on the TNF-α campaign.

### Consequences for the campaign

- The gate rejects, it does not select. A design scoring 0–1 features across samples can be
  discarded; a design scoring 3–4 has **not** been shown to be an agonist.
- Report mean across samples, never best-of-n: `max features` would have called exendin-4 the
  best ligand in the panel on a single lucky sample.
- Short binders need a separate sanity check, since the one false positive here was the shortest
  ligand tested.
- An independent predictor is required before any ranking is believed. Chai-1 is installed and
  its score capture is working.
