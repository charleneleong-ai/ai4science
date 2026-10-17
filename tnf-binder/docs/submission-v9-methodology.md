# Conditional TNF-α binder — methodology

Anthropic × Adaptyv Challenge 02, Track 3. Twenty de novo single-chain designs.

## About

| | |
| --- | --- |
| **Designs submitted** | 20 |
| **Experimentally validated** | 0 — nothing here is measured binding |
| **Design class** | minibinder 100% (40–100 aa); 18 × 91 aa, 2 × 74 aa |
| **Classification** | mainly alpha 90% (helix_0 chassis), mainly beta 10% (free_1 chassis) |
| **De novo** | 100% — generated from scratch, no binder used as a starting point |
| **Target** | human TNF-α soluble homotrimer, UniProt P01375-1, residues 77–233 |
| **Epitope** | conserved A/91 lobe (A86–92 + B33/B147), 90% species-conserved |

## Objectives, and how they conflict

The challenge ranks pH-selectivity first, mouse cross-reactivity second, affinity third. The
conflict between them is concentrated in one residue.

The epitope was computed rather than taken from the literature: TNF residues within 4.5 Å of
TNFR2 in **3ALQ**, intersected with a human/mouse alignment of the soluble domain (124/157 =
79.0% identity). That gives **33 contact positions, 24 conserved and 9 divergent**.

**R31→Q is the trap.** It carries the largest buried area in the entire epitope (193 Å²) and sits
beside R32 in the 29–33 loop that adalimumab and infliximab engage — which is why neither drug is
mouse cross-reactive. Any design gripping that loop wins objective 1 and loses objective 2.

Clustering the conserved exposed contacts into compact patches gave a better target:

| epitope centre | conserved / divergent SASA | conserved | pH anchor |
| --- | --- | --- | --- |
| **A/91** — A86-92 + B33/B147 | 685 / 72 Å² | **90%** | **K90** |
| B/147 — B20-33, B143-149 | 573 / 68 Å² | 89% | R32 |
| A/115 — 63, 67, 113, 115, 143-149 | 378 / 0 Å² | 100% | none |

The **A/91 lobe** was chosen: a genuine functional receptor site, 90% species-conserved by
construction, carrying conserved **K90** as a cationic anchor for the pH switch. One epitope
therefore serves all three objectives. Chain C sits 18.7 Å away, so the target was truncated to
the A+B protomer pair (304 residues rather than 471).

## The pH switch, and three hypotheses that failed

The intuitive mechanism — place a histidine beside K90 so it protonates at pH 6.0 and repels — is
self-defeating. A neighbouring cation **lowers** histidine pKa, often to 5.0–5.5, so the histidine
stays neutral at pH 6.0 and the switch never fires. Pairing histidine with Asp/Glu raises pKa but
makes protonation *stabilising*, which is the wrong direction entirely.

Three mechanistic hypotheses — buried-histidine desolvation, cation repulsion, carboxylate
proximity — were each tested against the pKa data and **each refuted**. No validated mechanism
replaced them, so positions come from a brute-force scan plus a physics solver (PROPKA on the
mutant in its bound context) rather than from a rule.

Every combination was measured **as a combination** rather than composed from singles, which was
necessary: `A35H+A38H` drops A35 from 6.07 to 5.94, out of the window, and `A70H+A72H` drops A70
from 6.44 to 6.23.

## The main methodological result: sampling depth

The pre-registered criterion was "predicted pKa in the 6.0–7.6 window on **every** sample",
measured on **three** diffusion samples. Re-measuring 25 designs at **twenty** samples — 500
structures, one PROPKA run per histidine per structure — changed the answer three ways.

### Two of the four positions are not switches

| His | designs | mean in-window rate | median pKa spread |
| --- | --- | --- | --- |
| **67** | 13 | **0.94** | 1.56 |
| 72 | 11 | 0.84 | 1.98 |
| 70 | 12 | 0.83 | 0.45 |
| 21 | 13 | **0.56** | **4.28** |
| 38 | 7 | 0.56 | 1.89 |
| 35 | 6 | **0.43** | **5.70** |

A median spread of 4.28 and 5.70 pH units across structures the folding model scores as equally
good is not a measurement with a window in it — the predicted pKa is a property of the sample,
not of the design.

### The criterion itself had to change

"In window on every sample" is a **maximum**, so it gets strictly harder as samples are added and
cannot be compared across sample counts. At n=3 most designs passed a ≥2-histidine gate; at n=20,
**0 of 25** do. The comparable quantity is the **in-window rate**, which does not depend on n.

The 0.90 rate threshold is **post hoc**. The full sensitivity curve is reported rather than a
single cut:

| in-window rate ≥ | designs with ≥2 such His | with ≥1 |
| --- | --- | --- |
| 1.00 (all samples) | 0 | 7 |
| 0.95 | 2 | 16 |
| **0.90** | **5** | 18 |
| 0.85 | 7 | 19 |
| 0.80 | 9 | 19 |
| 0.75 | 10 | 22 |

### The rate is design-specific, not positional

His70 and His72 hold up individually but degrade to **0.65** and **0.60** inside
`E67H_E70H_R72H`. Switch counts therefore cannot be composed from single-substitution scans.

### Mouse retention moved in both directions

Retention moved by up to **0.25 in both directions** between 3 and 20 samples, which reversed
which designs ranked best. Because the submission order turned on it, the top-tier designs still
carrying 3-sample mouse data were re-folded at 20 samples. All three 3-sample values were
optimistic — 0.874 → **0.802**, 0.824 → **0.735**, 0.853 → **0.860** — and the order changed.
All five top-ranked rows now rest on 20 samples.

## Ranking method

Rows are ordered by the challenge's own objective order, each term at the deepest sampling
available:

1. **Switching histidines** — a histidine counts if its predicted pKa sits in the 6.0–7.6 window
   on ≥90% of 20 samples. Within an equal count, a design measured at 20 samples is placed above
   one still provisional on 3, so an unmeasured design is neither credited nor penalised as
   though it had been tested.
2. **Mouse retention** — the ratio of symmetrised binder-to-target ipTM means, mouse over human.
3. **Human binder-to-target ipTM** — the **mean** over samples, not best-of-n.

| # | name | design | switches | mouse retention | samples |
| --- | --- | --- | --- | --- | --- |
| 1 | tnf2_15 | `R35H_E67H_R72H` | 2 | **1.041** | 20 |
| 2 | tnf2_13 | `R35H_E67H_E70H` | 2 | 0.893 | 20 |
| 3 | tnf2_05 | `E67H_E70H` | 2 | 0.860 | 20 |
| 4 | tnf2_03 | `E67H_R72H` | 2 | 0.802 | 20 |
| 5 | tnf2_01 | `R21H_E67H_E70H` | 2 | 0.735 | 20 |
| 6 | tnf2_12 | `R21H_R35H_E67H` | 1 | 1.085 | 20 |
| 7–16 | | helix_0 series | 1 | 0.726–0.951 | 3–20 |
| 17–18 | tnf2_06, tnf2_07 | | 1, provisional | 0.689–0.741 | 3 |
| 19–20 | tnf2_10, tnf2_09 | free_1 chassis | 0 | 0.482, 0.466 | 20 |

## Measurement discipline

Interface confidence is **not** reported as a point estimate — on identical structures AF2 scored
a BindCraft design **0.92** where Boltz-2 scored **0.21**.

Two artefacts were found and corrected, and both would have inverted conclusions:

- **Global ipTM is unusable on a multi-chain target.** The TNF trimer's own internal interfaces
  score ~0.88 and swamp the binder's. Moving dimer → trimer lifted global ipTM 0.20 → 0.52 while
  the binder-to-target interface stayed flat. Every number here is **binder-to-target only**.
- **Best-of-n inflates.** Metrics are stored per diffusion sample; ranking uses the **mean**.

No model scores designs it also optimised. And one design was nearly discarded on a single model
reading **232 Å²** with a tiny-interface flag, against a 20-sample median of **825 Å²** that never
fell below 535 — which is why single-model geometry is reported alongside the 20-sample
distribution rather than instead of it.

## Funnel

Failed generations are in the denominator.

| stage | n |
| --- | --- |
| BoltzGen designs generated | 100 |
| campaign constructs (parents + histidine variants) | 57 |
| folded parents | 16 |
| chassis passing **both** the fold-novelty and composition gates | **2** |
| designs folded at 20 samples | 25 human / 28 mouse (500 + 560 structures) |
| designs with ≥2 switching histidines at rate ≥0.90 | **5** |
| designs with ≥1 | 18 |
| submitted | 20 |

The concentration on two chassis is a **result of the funnel, not a preference**: of 16 folded
parents, only `screen_pocket_helix_0` and `screen_pocket_free_1` passed both gates. 18 of the 20
rows are a measured histidine series on the first chassis, pairwise identity **0.934–0.989**, so
they are a substitution series rather than independent designs.

## Stack

| role | tool | version / checkpoint |
| --- | --- | --- |
| design (sampling) | BoltzGen | 0.1.4, `boltzgen1_adherence` + `boltzgen1_diverse` |
| design (optimisation) | BindCraft | 1.5.3, AF2 `params_2022-12-06` |
| design (dual-target) | BindCraft2 | 1.0.3, jax 0.11.2 |
| co-folding / validation | Boltz-2 | 2.2.1, `boltz2_conf.ckpt` |
| side-chain building | PyRosetta | 2026.29 |
| pKa prediction | PROPKA | 3.5.1 |
| sequence novelty | MMseqs2 | vs Swiss-Prot + PDB seqres, snapshot 2026-10-10 |
| fold novelty | FoldSeek + TM-align | CATH50, snapshot 2026-10-11 |

Boltz-2's affinity head is **not used anywhere**: it requires the binder chain to be a
small-molecule ligand and caps at 128 atoms (trained at ≤56), against ~450 heavy atoms for these
designs.

## What is not established

- **No measured binding of any kind.** No expression, no SEC, no thermal stability, no SPR or BLI.
- pH selectivity rests on PROPKA over **predicted** structures. No folding model represents
  protonation, and the pH 7.4 / 6.0 assay is not simulated.
- Mouse cross-reactivity rests on a **predicted** mouse trimer and conserved-epitope coverage, not
  on an experimental mouse structure.
- Three mechanistic hypotheses for the switch were each refuted by the pKa data; positions come
  from a scan, not a rule, and no validated mechanism replaces them.
- **The filter set was not fixed before results were seen.** The composition gate, the fold gate
  and a lowered burial floor (15 Å² → 8 Å²) were each added in response to an observed failure.
  The funnel above is therefore not a clean pre-registered statistic. Every deviation is declared
  in the experiment card.
- Two designs (`tnf2_06`, `tnf2_07`) were never folded at 20 samples; their switch counts are
  provisional on 3 and they are ranked below every design measured at 20 with the same count.

## Data

Full code, per-design metrics, every structure and the complete pKa tables:
<https://github.com/charleneleong-ai/ai4science/tree/main/tnf-binder>
