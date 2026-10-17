# Conditional TNF-α binder

20 de novo single-chain minibinders. No binding is measured. The headline is a negative result
about our own pH-switch evidence, found by verifying it three ways.

## Epitope, computed not cited

TNF residues within 4.5 Å of TNFR2 in 3ALQ, intersected with a human/mouse alignment: 33
contacts, 24 conserved, 9 divergent. **R31→Q is the trap** — largest buried area in the epitope
(193 Å²), in the loop adalimumab and infliximab engage, which is why neither is mouse
cross-reactive. Any design gripping it wins objective 1 and loses objective 2. Clustering
conserved exposed contacts gave the **A/91 lobe**, with conserved **K90** as a cationic anchor.
Three switch hypotheses — buried-His desolvation, cation repulsion, carboxylate proximity — were
each tested and **each refuted**, so positions come from a scan, not a rule.

## Three verifications, which disagree

**1. Depth.** The pre-registered criterion was pKa in-window on every sample, at three samples.
At twenty (500 structures) His21 and His35 showed median pKa spreads of **4.28** and **5.70** pH
units between structures the model scores equally. "Every sample" is also a *maximum*, so it
tightens as n grows and cannot be compared across sample counts; converted to a rate. The 0.90
cut is post hoc, so the curve is reported: ≥1.00 → 0 designs with two switches, 0.95 → 2,
**0.90 → 5**, 0.85 → 7, 0.80 → 9.

**2. An independent method.** Every claim rested on PROPKA, an empirical scheme; re-running it
on more structures measures structural variability, not method bias. **pKAI** (graph neural net)
shares no machinery with it and reads systematically **+0.76 higher**. The two disagree on
**23 of 62** design-positions at the 0.90 cut.

**3. Coupling to binding — this inverted the result.** A histidine whose pKa is the same bound
and free contributes *nothing* to pH-dependent binding: the protonation equilibrium cancels
between states. What governs pH-dependence of Kd is the linkage, Δ = pKa(complex) − pKa(free),
which must be negative for protonation to drive dissociation. Never previously measured.

| His | in-window rate | mean Δ | coupled samples |
| --- | --- | --- | --- |
| 67 | **0.94** | −0.09 | **0/13** |
| 70 | 0.83 | −0.05 | **0/12** |
| 72 | 0.84 | −0.33 | 6/11 |
| 21 | 0.56 | **−0.49** | **12/13** |
| 35 | 0.43 | **−0.82** | **6/6** |
| 38 | 0.56 | −0.39 | 3/7 |

**The two analyses point opposite ways.** His67 and His70 — the best-determined positions, which
every top design was built on — are uncoupled on all 25 samples: they titrate cleanly inside the
window and that does nothing to binding. His21 and His35, which depth *rejected*, are the most
strongly coupled.

Geometry does not explain it away: His67 (9.7 Å from target) and His21 (9.0 Å) are equidistant,
yet one is uncoupled on every sample and the other coupled on 12 of 13. But the distances matter
— **all six histidines sit 9–21 Å from the target, none in direct contact** — the likeliest
reason coupling is weak throughout. **Only His72 passes all three checks, in two designs.** Δ is
rigid-body: the free state is the binder chain taken from the complex, not relaxed.

## Ranking

Graded tiers, so a check that rejects everything degrades the order rather than flattening it.
**A** = both pKa methods in-window on ≥90% of 20 samples *and* coupled; **B** = both agree,
coupling too weak to call; **C** = one method only, reported but **not** promoting. Then mouse
retention, then human ipTM as a mean, not best-of-n. Retention itself moved up to 0.25 in both
directions between 3 and 20 samples, so the top rows were re-folded at 20 before ordering.
`tnf2_15` and `tnf2_03` hold the only tier-A switches, both His72.

## Funnel

```
100 generated -> 57 constructs -> 16 folded parents -> 2 chassis passed both
fold-novelty and composition gates -> 25 folded at 20 samples (1060 structures)
-> 5 carry >=2 switching His by PROPKA alone -> 2 survive all three checks
-> 20 submitted
```

Failed generations are in the denominator. 18 of 20 rows are a histidine series on one chassis
(identity 0.934–0.989) — a funnel result, not a preference.

## Liabilities

**18 of 20 designs carry a single unpaired cysteine at position 34**, from the helix_0 chassis.
A free thiol risks disulfide scrambling, covalent dimerisation and aggregation. It was not
designed out and no cysteine-free variant is verified. On 2,214 wet-lab-tested public designs the
best expression predictor we found was glutamate fraction (AUROC 0.724): helix_0 at 0.283 is
favourable, free_1 at 0.014 is not.

## Not established

pH selectivity rests on predicted pKa over predicted structures; the assay is not simulated.
Mouse cross-reactivity rests on a predicted mouse trimer. Filters were **not** fixed before
results were seen. On the evidence above the pH-switch claim for most of these designs is
**weaker than a single-method analysis would suggest**, and we report that.

github.com/charleneleong-ai/ai4science/tree/main/tnf-binder
