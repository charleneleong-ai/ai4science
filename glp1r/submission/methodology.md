# De novo GLP1R binders

20 single-chain minibinders, 43–91 aa, against the extracellular-domain groove. Nothing is
measured. **We do not claim these are agonists** — the evidence says the opposite, and that is
the most useful thing we can report.

## The target, measured not cited

Receptor residues within 4.5 Å of GLP-1 in **6X18**: 52 contacts, splitting cleanly into **29
for the peptide's C-terminal half** (19 in the ECD) and **26 for its N-terminal half** (all in the
TM bundle and loops). Only three see both. The site centroids are **31.6 Å apart** and the peptide
N-terminus sits **19.2 Å below** the vestibule rim, so a molecule that binds *and* activates must
span the receptor and terminate deep inside the bundle. The folded ECD is ordered only 29–128
across all 13 structures resolving it, so that is the design target.

## Why we expect these to be antagonists

Six independent campaigns have put a protein in this groove, including a 10¹⁰ antibody library
carrying >150,000 GPCR-binding motifs in HCDR3. **Zero full agonists**; the best ECD-restricted
agonism reported is Emax <2% of exendin-4. A 10,080-design de novo campaign against this exact
epitope, read out in a CRE-luciferase line, produced two antagonists at 39 and 61 nM.

The activating pharmacophore is GLP-1 7–17 inside the bundle: an 11-mer with no ECD-binding helix
at all is a 0.5 nM full agonist. The groove is the affinity trap, not the switch.

## Verification: three checks, declared before generating

**1. The activation gate was calibrated on known-function ligands, and it failed.** A gate of
six pharmacophore features plus insertion depth orders *experimental* structures correctly —
GLP-1 6/6, exendin-4 and exendin-P5 4/6, the weak partial GLP-1(9-36) 1/6, a groove-only
truncation 0/6. On **Boltz-2 predictions** of the same ligands it collapses: an **inactive 15-mer
ranks second**, above two full agonists, because with no N-terminal half nothing anchors it to
the groove and the model drops it into the bundle. The gate **rejects but does not select**, and
short designs are at specific risk. One sample gave native GLP-1 0/6 against its experimental
6/6; only all five recovered signal.

**2. Nothing is ranked by the model that made it.** BoltzGen generates and ranks with Boltz-2
weights. Every survivor was re-folded with **Chai-1**, which shares no weights with it.

**3. Scores were calibrated against wet-lab outcomes.** On 2,214 publicly measured designs,
within-target AUROC for binding: interface pLDDT 0.652, ipTM 0.607 — nothing reaches 0.66, and
every score predicts *expression* better than binding.

## Funnel

```
200 designs, 4 hotspot variants on the ECD core (29-128)
140 pass sequence triage (70%)   <- 46 alanine over-representation, 27 Ala+Gly
 23 also reach BoltzGen iPTM >=0.5
 21 confirmed by Chai-1 at >=0.65
 20 submitted (cap)
```

All four variants cleared the 0.11–0.30 interface plateau a prior campaign hit on a buried
pocket; this groove is an easier target.

## Two models agree they bind; neither supports ranking them

| | BoltzGen | Chai-1 |
| --- | --- | --- |
| median | 0.535 | **0.802** |
| range | 0.502–0.616 | 0.150–0.862 |

Spearman between them is **+0.377** — the 21 survivors span only 0.651–0.862 on the independent
model, so they are not distinguishable and **no fine-grained ranking is claimed**. Row order is by
Chai-1, the model that did not generate them, and should be read as weak.

The disagreement did real work: two designs scoring 0.551 and 0.516 on BoltzGen fell to **0.236**
and **0.150** on Chai-1 and were dropped. Both would have made a top-20 ranked on the generator
alone.

Mutual identity among the 20 is **0.011–0.317, median 0.096** — a diverse set, not a
substitution series.

## Novelty, checked

Sequence: MMseqs2 against SwissProt + PDB seqres (1.75M sequences). Every one of the 20 clears
the 30% identity condition on any significant measure — **the best hit for any design has
e-value 0.27, and most exceed 10**, so these are alignment noise rather than homology.

Structure: Foldseek against CATH50, scoring the **design chain only**. **19 of 20 fall in the
"moderate" band** (TM 0.54–0.75 at >70% coverage), which a prior campaign saw pass this scorer.
**One — the last row — is HIGH (TM 0.844 at full coverage)** and may be rejected by the de novo
filter. It is ranked last deliberately so it is the row dropped; its best alternative scored
worse (TM 0.869), so no swap was available.

## Not established

No measured binding, expression, or activation. Two predictors agreeing raises confidence but is
not a measurement — on the public corpus above, ipTM's own AUROC for binding is 0.607. **Agonism
is explicitly not claimed**: these engage the ECD groove, the published route to an antagonist,
and our own activation gate cannot select for agonism on predicted structures. The novelty
thresholds (≤30% identity, sub-moderate structural similarity) have not been independently
verified against the organisers' scorer.
