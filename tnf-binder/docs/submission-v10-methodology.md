# Conditional TNF-α binder

20 de novo single-chain minibinders. No binding is measured. The headline is a negative result
about our own pH-switch evidence, found by verifying it three ways.

## Epitope, computed not cited

TNF residues within 4.5 Å of TNFR2 in 3ALQ, intersected with a human/mouse alignment
(124/157 = 79.0% identity): 33 contacts, 24 conserved, 9 divergent. **R31→Q is the trap** —
largest buried area in the epitope (193 Å²), in the 29–33 loop adalimumab and infliximab engage,
which is why neither is mouse cross-reactive. Any design gripping it wins objective 1 and loses
objective 2. Clustering conserved exposed contacts gave the **A/91 lobe**, with conserved **K90**
as a cationic anchor. Three mechanistic hypotheses for the switch — buried-His desolvation,
cation repulsion, carboxylate proximity — were each tested and **each refuted**, so positions
come from a brute-force scan, not a rule.

## Three verifications, which disagree

**1. Depth.** The pre-registered criterion was pKa in the 6.0–7.6 window on every sample, at
three samples. At twenty (500 structures) His21 and His35 showed median pKa spreads of **4.28**
and **5.70** pH units between structures the model scores equally. "On every sample" is also a
*maximum*, so it tightens as n grows and cannot be compared across sample counts; converted to a
rate. The 0.90 cut is post hoc, so the curve is reported: ≥1.00 → 0 designs with two switches,
0.95 → 2, **0.90 → 5**, 0.85 → 7, 0.80 → 9.

**2. An independent method.** Every claim rested on PROPKA, an empirical scheme with fitted
coefficients; re-running it on more structures measures structural variability, not method bias.
**pKAI** (graph neural net) shares no machinery with it and reads systematically **+0.76
higher**. The two disagree on **23 of 62** design-positions at the 0.90 cut.

**3. Coupling to binding — the one that inverted the result.** A histidine whose pKa is the same
bound and free contributes *nothing* to pH-dependent binding: the protonation equilibrium cancels
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
every top design was built on — are uncoupled on all 25 samples. They titrate cleanly inside the
window and that titration does nothing to binding. His21 and His35, which depth *rejected*, are
the most strongly coupled.

Geometry does not explain it away: His67 (9.7 Å from target) and His21 (9.0 Å) are equidistant,
yet one is uncoupled on every sample and the other coupled on 12 of 13. But the distances matter
— **all six histidines sit 9–21 Å from the target, none in direct contact** — which is the
likeliest reason coupling is weak throughout. **Only His72 passes all three checks, in two
designs.** Δ is rigid-body: the free state is the binder chain extracted from the complex, not
relaxed, so it isolates the target's electrostatic contribution only.

## Ranking

Graded tiers, so a check that rejects everything degrades the order rather than flattening it.
**A** = both pKa methods in-window on ≥90% of 20 samples *and* coupled; **B** = both agree,
coupling too weak to call; **C** = one method only, reported but **not** promoting — a switch
only one method believes in is the claim this exercise exists to stop making. Then mouse
retention, then human ipTM as a mean, not best-of-n. Retention itself moved up to 0.25 in both
directions between 3 and 20 samples, so the top rows were re-folded at 20 before ordering.

Result: `tnf2_15` and `tnf2_03` hold the only tier-A switches, both His72.

## Funnel

```
100  BoltzGen designs generated
 57  campaign constructs
 16  folded parents
  2  chassis passed BOTH fold-novelty and composition gates
 25  folded at 20 samples (500 human + 560 mouse structures)
  5  carry >=2 switching His by PROPKA alone
  2  carry a switch surviving all three verifications
 20  submitted
```

Failed generations are in the denominator. 18 of 20 rows are a histidine series on one chassis
(identity 0.934–0.989) — a funnel result, not a preference.

## Not established

No measured binding. pH selectivity rests on predicted pKa over predicted structures; the assay
is not simulated. Mouse cross-reactivity rests on a predicted mouse trimer. Filters were **not**
fixed before results were seen. On the evidence above the pH-switch claim for most of these
designs is **weaker than a single-method analysis would suggest**, and we report that rather
than the single-method number.

<https://github.com/charleneleong-ai/ai4science/tree/main/tnf-binder>
