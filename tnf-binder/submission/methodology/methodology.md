# Conditional TNF-α binder — methodology

## This is a designed series, not 12 independent guesses

The set is built to be **informative whichever way it reads in the lab**. Ten designs are a
histidine titration on one backbone spanning predicted pKa 6.07–7.12, with one, two or three
histidines; two more are a **matched pair on a second backbone — `tnf2_01` carries the switch,
`tnf2_12` is the identical sequence with zero histidines.** That pair is the control: if
`tnf2_01` shows pH-dependent release and `tnf2_12` does not, the switch is the cause, and if
neither does, pKa-engineered switching is ruled out on this scaffold rather than left ambiguous.

## Method: what is different here

Everything below is open — code, tests, the full funnel including every failure, and the
measurements behind each number:
**https://github.com/charleneleong-ai/ai4science/tree/feat/tnf-submission-2-results/tnf-binder**

The per-design metrics are [`docs/submission-v2-evidence.csv`](https://github.com/charleneleong-ai/ai4science/blob/feat/tnf-submission-2-results/tnf-binder/docs/submission-v2-evidence.csv),
the raw pKa-across-structures measurements are [`docs/pka_robustness.csv`](https://github.com/charleneleong-ai/ai4science/blob/feat/tnf-submission-2-results/tnf-binder/docs/pka_robustness.csv),
and the gates are [`src/tnf_binder/verifiers/`](https://github.com/charleneleong-ai/ai4science/tree/feat/tnf-submission-2-results/tnf-binder/src/tnf_binder/verifiers)
with 93 tests, each anchored to an observed value rather than a synthetic one.

**1 · A verifier stack ported from a different design problem.** The gates reuse the `Verdict`
and `Verifier` contract of [`touchstone`](https://github.com/charleneleong-ai/ai4science/tree/main/touchstone),
a verifier we built for *metal-binding* proteins. Its metal geometry does not transfer, but its
architecture does — in particular the distinction between **weak** (judged and failed) and
**defer** (cannot judge). Conflating the two let a missing measurement read as a pass; keeping
them separate is what stopped a design with no pKa prediction being scored as if it had one.

**2 · Tiers drop candidates; they never down-rank them.** The top-scoring construct of an
earlier round was verbatim human ubiquitin, and the top of a 40-design shortlist was 59 %
alanine. No ranking function rescues either — only a gate does.

**3 · A switch-robustness tier.** Every pKa was re-measured against all three diffusion samples
of its own complex, not just the structure it was designed on. **8 of 23 histidines do not hold
their window across all three.** One reads 3.99–6.70 depending on which structure you ask. We
have not seen this check applied to designed pH switches, and it reordered our entire ranking.

**4 · Ranking on a physical quantity, not a feature count.** ΔQ — the interface charge gained
between the two assay pH values — evaluated at each histidine's *least favourable* structure.
A lower bound, not a best case.

**5 · Combinations measured, never composed.** All 41 histidine combinations were scanned,
because single-mutant pKa values do not add: `A35H+A38H` pushes A35 out of the window.

**6 · The negative results are in the repository.** Eight hypotheses we held were refuted by
measurement, including three separate mechanisms for *why* a histidine should switch — burial,
cation repulsion, carboxylate proximity. All three have the wrong sign or no effect. The
positions used here are empirical, from an exhaustive scan plus a physics solver, because every
mechanism we reasoned our way to was wrong.

**7 · One claim was retracted.** We reported that our fold-novelty gate reproduced Proteinbase's
verdicts 3/3. It did not: two of the three only scored correctly because a missing identity
value defaulted to zero. The gate is conservative, not calibrated — it under-scores, which your
own check confirmed when it returned 4/4 for a design we scored "moderate". The retraction and
a regression test pinning the bug are in the repository.

## Epitope: computed, not taken from literature

TNF-α residues within 4.5 Å of TNFR2 in **3ALQ**, intersected with a human/mouse alignment of
the soluble domains (**124 of 156** aligned residues, 79 %), give **33 contact positions: 24 conserved,
9 divergent** (P20H, R31Q, S71D, T72–, H73Y, V85I, T89E, I97V, A111L).

### Cross-reactivity was then measured, not inferred

Conserved-epitope coverage is a proxy, so each design was additionally **folded against the
mouse TNF-α trimer** (P06804, the sequence on the challenge page) and its interface compared
with the human one under the identical protocol — a within-model ratio, which is what this
comparison can legitimately support.

**9 of 12 retain ≥70 % of their human interface confidence on mouse**, median 86 %, and four
score higher on mouse than on human. The proxy and the direct prediction disagree in one place
that matters: the design with among the highest conserved-epitope coverage (0.79) retains only
**47 %** on mouse, and its control 41 %. Both were demoted. Designing against the conserved
patch makes cross-reactivity likely; it does not make it true, and only folding against the
other ortholog distinguishes the two.

**R31→Q is the trap.** It carries the largest solvent-accessible area in the epitope (193 Å²)
and sits beside R32 in the 29–33 loop that adalimumab and infliximab both engage — which is
why neither is mouse cross-reactive. Those nine positions were set as coldspots and the
conserved patches targeted instead, so cross-reactivity is designed in rather than hoped for.

## pH switch: scanned, then measured in combination

PROPKA 3.5.1 over PyRosetta-built side chains, on the **bound complex**, not the free binder.
Every interface position was substituted to histidine in turn (23 on the lead chassis), giving
6 switchable sites, 4 of them in the ideal 6.3–7.2 window.

**Combinations were measured, never composed from single-mutant values.** This mattered:
`A35H+A38H` drops A35 from 6.07 to **5.94, out of the window**, and `A70H+A72H` drops A70 from
6.44 to 6.23. All 41 combinations of the six sites were scanned.

**And every surviving pKa was then re-measured against three independent predicted structures**
rather than the one it was designed on. This is the check we would most want someone else to
have run on us: 8 of 23 histidines do not hold their window across all three, and `R21H` alone
ranges 3.99–6.70 where its design structure read 7.05. Five designs hold every histidine in
window on every structure; they lead the ranking, and the single most reproducible switch is
`R72H` at 6.37 with a spread of 0.27.

Three mechanisms we initially reasoned toward were each refuted by the scan: burying the
histidine *lowers* its pKa (r = −0.37), a neighbouring cation has the wrong sign (r = +0.11),
and carboxylate proximity gives only r = +0.21 with counter-examples at net +6 sitting at
pKa 3.6. The positions here are empirical, from an exhaustive scan plus a physics solver.

## The verifier stack

Each tier judges something the tier above cannot see, and tiers **drop** candidates rather than
down-rank them. Counts are for this round: 16 folded parent designs in, 12 submitted.

```mermaid
flowchart TB
    G["16 folded parent designs<br/>BoltzGen screens"] --> V1
    subgraph S["verifier stack — each tier drops, never down-ranks"]
        direction TB
        V1["<b>1 · composition</b><br/>model-free<br/><i>no model can be trusted to<br/>reject its own degenerate output</i>"]
        V2["<b>2 · fold novelty</b><br/>FoldSeek + TM-align vs CATH50"]
        V3["<b>3 · interface</b><br/>buried area · clashes · epitope focus"]
        V4["<b>4 · cross-reactivity</b><br/>interface area off the 9 divergent positions"]
        V5["<b>5 · pH switch</b><br/>PROPKA on the bound complex"]
        V6["<b>6 · switch robustness</b><br/>same pKa across 3 predicted structures"]
        V7["<b>7 · switch magnitude</b><br/>ΔQ at the interface, pH 7.4 → 6.0"]
    end
    V1 -->|"9 survive"| V2
    V1 -.->|"7 killed · poly-Ala 50-75%"| X1[["eliminated"]]
    V2 -->|"2 chassis"| V3
    V2 -.->|"7 killed · high fold similarity"| X2[["eliminated"]]
    V3 -->|"15 His variants built"| V4
    V4 -->|"12 submitted"| V5
    V4 -.->|"3 killed · <0.6, incl. best interface<br/>at 59% on divergent residues"| X4[["eliminated"]]
    V5 -->|"41 combinations measured"| V6
    V6 -->|"5 hold every His<br/>in window on all 3"| V7
    V6 -.->|"8 of 23 His are<br/>structure-dependent"| DEMOTE[["demoted, not dropped"]]
    V7 -->|"ranked on the<br/>guaranteed swing"| OUT["ranked submission"]
```

**Tier 6 is the one that changed the answer.** Every pKa originally came from a single predicted
structure. Re-scoring each histidine against all three Boltz diffusion samples showed **8 of 23
are structure-dependent** — including the design this submission originally led with, whose
switch reads 6.27 on one structure and 5.92 averaged over three. The ranking below is the
corrected one; nothing was dropped on this basis, because a structure-dependent prediction is
weak evidence rather than evidence of failure.

## Selection

Designs passed independent gates before ranking, and the gates **dropped** candidates rather
than down-ranking them — the top-scoring construct of an earlier round was verbatim ubiquitin,
and the top of a 40-design shortlist was 59 % alanine.

- **composition** — model-free; rejects low-complexity and poly-alanine output
- **fold novelty** — FoldSeek + TM-align against CATH50
- **interface** — buried surface area and steric clashes, not contact counts
- **cross-reactivity** — fraction of buried interface area *off* the 9 divergent positions

### Ranking: the swing, not the count

Counting "histidines in the window" treats a pKa of 6.0 and 6.7 as equal when they are not, so
the ranking uses the quantity the assay actually reads: **ΔQ, the charge the interface gains
going from pH 7.4 to pH 6.0**, by Henderson–Hasselbalch on each measured pKa. It peaks at
pKa 6.7, the midpoint of the two assay points — so distance from 6.7, not a low pKa as such, is
what costs a design.

Each design is ranked on its **ΔQ floor**: every histidine scored on whichever of its three
structures gives the smaller swing. That is a lower bound on the switch, not a best case.

| | ΔQ floor | robust His |
| --- | --- | --- |
| `tnf2_01` (3 His) | **1.38** | 3/3 |
| `tnf2_04`, `tnf2_05` (2 His) | 0.89–0.99 | 2/2 |
| `tnf2_09` (1 His) | 0.57 | 1/1 — most reproducible single, spread 0.27 |
| `tnf2_12` (1 His) | **0.01** | 0/1 — the switch vanishes on its worst structure |

One ideally-placed histidine can contribute at most ~0.62, so the three-histidine lead design
is carrying roughly 2.2 histidines' worth of guaranteed swing.

Objective 3 asks for **no detectable binding at pH 6.0**, not merely weaker binding, and that
is a higher bar than a shifted affinity. On that reading a ΔQ near 0.4 is unlikely to abolish
binding outright, and only the multi-histidine designs at the top of the table carry enough
guaranteed swing to be plausible candidates for a full switch. The single-histidine designs are
included to locate the threshold, not because we expect each to abolish binding.

`tnf2_06` is ranked above its ΔQ (0.37) deliberately: it is the only design with a credible
measured interface (923 Å² buried, zero cysteines, novelty 4/4), and **a design that does not
bind cannot show pH-dependent release at all**. `tnf2_07` follows it as its matched control.

Cross-reactivity, not affinity, did the discriminating: it ranged 0.41–0.89 across variants. One
design with the family's best interface confidence put **59 % of that interface on
human/mouse-divergent residues** and was dropped, along with two others below 0.60.

A measured trade-off worth stating: **each added histidine costs interface confidence** — mean
binder-to-target ipTM 0.235 for singles, 0.180 for doubles, 0.158 for triples. The series spans
that trade-off deliberately rather than picking one point on it.

## What is not established

**No measured binding of any kind.** No expression, SEC, thermal stability, SPR or BLI; every
number above is computational. pH selectivity rests on predicted pKa, and no folding model
represents protonation state — the robustness tier above bounds how much that prediction moves
with the input structure, which is not the same as showing it is right. Mouse cross-reactivity is a *predicted* interface against the
mouse sequence — a direct prediction rather than a coverage proxy, and the two disagree for one
design, but still not measured binding to mouse TNF-α.

Interface confidence is **model-dependent and unresolved**, so a range is reported rather than a
number: AlphaFold2 rated its own optimised designs at i_pTM 0.92 where Boltz-2 rated the same
structures 0.21. No design here was scored by a model that also optimised it. Global ipTM is not
used at all — on a homotrimer target it is inflated by +0.405 on average by the target's own
internal interfaces; all values are binder-to-target chain-pair ipTM, symmetrised and averaged
across diffusion samples rather than taken best-of-n.

One earlier design failed your novelty check at 37.5 % sequence identity while two passed at
38.1 %: the failure was **fold**, not sequence — an immunoglobulin-like β-sandwich at TM 0.893.
That rule shaped this round's selection.
