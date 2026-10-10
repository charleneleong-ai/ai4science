# CASMI 2026 spike — findings and whether to commit

> **Corrected 2026-10-08.** The first version of this writeup concluded that the public
> leaderboard "measures an exact join" and that a train-structure pool was the winning design.
> Both were wrong. A submission built on that premise scored **0.115 (rank 2397/2686)**
> against a locally predicted 0.664. The corrections are recorded below rather than silently
> overwritten, because how the error survived is the most useful thing here.

**Verdict:** not a realistic medal push. The field sits at 0.42–0.48 by retrieving from
PubChem; this work retrieved from train and scored 0.115. Closing that gap is a data-engineering
task — attach a PubChem-scale pool as a Kaggle Dataset — and is exactly what the leading public
notebooks already do. The fragment ranker may add something on top, but that is untested.

[Competition](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra):
up to 25 ranked SMILES per molecule from LC-MS/MS, scored MRR@25 on the tautomer-canonical
InChIKey connectivity block. Code competition: notebook only, internet off. Entry closes
2026-12-07.

## The real result

| | publicScore | rank |
|---|---|---|
| fragments over a train-structure pool | **0.115** | 2397 / 2686 |
| fragments over train + [PubChem tier](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-pubchem-tier) (106M), 4000 nearest-mass per molecule | **0.049** | — |
| train + 300 most-documented PubChem connectivities ([popularity prior](https://www.kaggle.com/datasets/dmitriigluzdov/casmi26-pubchem-popularity-prior)), fragments + 0.25 × popularity | **0.135** | — |
| same shortlist and prior, [public FPNet](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-fpnet-full1) fingerprint score in place of fragments | **0.250** | — |
| same, averaging two FPNet checkpoints (full1 + FPNet A) | **0.258** | — |
| same, popularity weight 0.15 instead of 0.25 | 0.254 | — |
| same as 0.258, shortlist 1000 instead of 300 | **0.261** | — |
| same, never proposing the train structure behind a copied spectrum | **0.262** | — |
| same, [ICEBERG](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-iceberg) re-ordering same-formula groups in the top 60 (CPU, 269 / 400 molecules scored) | **0.286** | — |
| same on a T4: ICEBERG scores all 371 molecules it can cover | **0.292** | — |
| same, [GLACIER](https://www.kaggle.com/datasets/ahmedberatozer/casmi26-glacier) fused beside ICEBERG (0.5 / 0.5) | **0.320** | — |
| leaderboard top | 0.480 | 1 |
| [public notebook, "v4n Fusion + PubChem"](https://www.kaggle.com/code/huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421) | 0.421 | — |

**A bigger pool alone made it worse.** With PubChem attached, 297 / 400 molecules have more
than 4000 structures inside 5 ppm; the cut keeps the nearest by mass, which is arbitrary among
same-formula isomers. Fragment explainability — 0.70 among ~46 train candidates — does not pick
the answer out of thousands. Run: 7767 s, PubChem mass convention matches RDKit exact mass
(median error 0.000000 Da). Which of the cut or the ranker loses more is not separable from one
score.

## What went wrong

**The local scorer was never calibrated against ground truth.** It scored submissions against
labels recovered by joining test spectra to bit-identical train spectra. Checked against the
0.421 public notebook only *after* submitting:

| submission | local scorer | real score |
|---|---|---|
| public notebook | 0.030 | **0.421** |
| this work | 0.664 | **0.115** |

A scorer that rates a known-good submission at 0.030 is not measuring the competition metric.

**The labels were mostly wrong.** Our top-1 equalled the recovered label for 215 of 400
molecules; had those been correct we would have scored at least 0.538. A 0.115 bounds them
correct for **at most ~46**.

**And the loop was closed.** Labels came from train, the pool *was* train, so the label was in
the pool by construction. "Recall 1.0000" was guaranteed, not measured — it was reported as a
ceiling.

## What is still true

- **The spectra really are duplicated.** 1189 / 1213 test spectra are bit-identical to a train
  row in both m/z and intensities. What is retracted is the inference drawn from it: the train
  structures attached to those spectra are mostly *not* the answers. **Why** remains unexplained.
- **A large share of answers lies outside train.** The 0.421 notebook's top-1 is present in
  train for only 12 / 400 molecules. This doesn't give the exact in-train fraction — that
  notebook is itself only partly right — but it rules out a train-only pool as a viable design.
- **The identity chain is sound.** Train's `inchikey14` equals RDKit's InChIKey14 from
  `normalized_smiles` on 3000 / 3000 rows checked. The failure was in labels, not keys.

## A valid result, correctly scoped

[`split-eval`](../src/casmi/cli.py) holds out *train* molecules and scores them against their
own structures, so its labels are genuine:

| ranker | MRR@25 |
|---|---|
| mass-error floor | 0.4119 |
| fragment explainability | **0.7011** |

That is a real ranking improvement — among ~46 candidates of which ~29 share the truth's
formula, with no training and no model weights. But it measures a regime where the answer is
guaranteed to be in the pool, which the competition does not offer. It does not transfer.

## Retracted claims

| claimed | status |
|---|---|
| the public LB measures an exact join; teams are "optimising a lookup" | **wrong** — the joined structures are mostly not answers |
| local leak-free estimate MRR@25 0.664 | **wrong** — scored against wrong labels in a closed loop |
| a train-structure pool gives recall 1.00 on the test set | **circular** — guaranteed by construction |
| COCONUT covers 0 / 400 test molecules | **no valid basis** — measured against the wrong labels |
| the test set is synthetic drug-like, 36% halogen | **no valid basis** — described the wrong labels; true answers' chemistry unknown |
| the duplication is a "data-prep error" | **unsupported** — duplication is a fact; intent was never knowable |

## Premises overturned along the way

| believed | measured |
|---|---|
| retrieval + rerank is the research frontier | de-novo generation is; graph diffusion first broke 0% exact-match on MassSpecGym |
| fragment enumeration cannot separate branching isomers | it can — a hydrogen-accounting bug in RDKit's fragmenter |
| the leaderboard's tight 0.43–0.48 band hid a field missing a free win | it was the honest difficulty of the task |

The last one is the most instructive. 2,685 teams clustered at 0.43–0.48 while I believed a
trivial join scored ~1.0. That contradiction was explained away instead of investigated.

## If continuing, in order

1. **Calibrate the scorer first.** Reproduce 0.421 for the public notebook locally before any
   other number is trusted. Labels must come from somewhere other than train.
2. ~~**PubChem-scale pool**~~ — done; 0.049 on its own.
3. ~~**A popularity prior**~~ — 0.135.
4. ~~**A learned spectral model**~~ — FPNet in place of fragments: 0.250, the largest single gain.
5. ~~**An FPNet ensemble**~~ — 0.258; a second checkpoint adds little.
6. ~~**Tuning the prior**~~ — weight 0.15: 0.254; shortlist 1000: 0.261. Both moves are
   within ~0.005 of 0.258, so the shortlist and the prior are no longer the bottleneck.
7. ~~**Excluding copied structures**~~ — 0.262. The 0.261 submission's top-1 was the train
   structure behind a copied spectrum for 267 / 400 molecules (FPNet full1 trained on those rows);
   the 0.421 notebook's was for 12 / 400, at median Tanimoto 0.05 to it. Excluding them changed
   all 267 top-1s and the score by 0.001: what replaced them is as rarely right. The copied
   structures also rule out library and analog search over train as built — its top hit is that
   structure, and the 0.421 answers are not its relatives.
8. ~~**A forward model**~~ — ICEBERG: 0.286, the second-largest gain after FPNet, and only
   269 / 400 molecules were scored (CPU, 1.6 predictions/s, 6 h budget; the GPU quota was spent).
   Changed 72 top-1s. Not separable from the switch to the pinned Python 3.12 image and RDKit
   2026.03 made in the same run: v11's output was not kept to compare against.

## What is reusable

The **mechanics**: the [notebook](../notebooks/kaggle_submission.py) runs in Kaggle's sandbox
in 266 s with rdkit installed offline from an attached wheel, and reproduces the package
pipeline byte-for-byte on all 400 rows; the [fragment ranker](../src/casmi/fragments.py) is
tested against isomers and ring cleavage.

The **lesson** is the larger asset. Every component here was verified — 93 tests, parity checks,
negative controls, bit-identity checks — and the composite was never checked against the one
external ground truth available for free: a public submission with a known score. That check
took two minutes and was done last.
