# Conditional TNF-α binder — methodology

## This is a designed series, not 12 independent guesses

The set is built to be **informative whichever way it reads in the lab**. Ten designs are a
histidine titration on one backbone spanning predicted pKa 6.07–7.12, with one, two or three
histidines; two more are a **matched pair on a second backbone — `tnf2_01` carries the switch,
`tnf2_12` is the identical sequence with zero histidines.** That pair is the control: if
`tnf2_01` shows pH-dependent release and `tnf2_12` does not, the switch is the cause, and if
neither does, pKa-engineered switching is ruled out on this scaffold rather than left ambiguous.

## Epitope: computed, not taken from literature

TNF-α residues within 4.5 Å of TNFR2 in **3ALQ**, intersected with a human/mouse alignment of
the soluble domain (124/157 = **79.0 %** identity), give **33 contact positions: 24 conserved,
9 divergent** (P20H, R31Q, S71D, T72–, H73Y, V85I, T89E, I97V, A111L).

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
6.44 to 6.23. All 41 combinations of the six sites were scanned; two triples hold 3/3
histidines in the ideal window.

Three mechanisms we initially reasoned toward were each refuted by the scan: burying the
histidine *lowers* its pKa (r = −0.37), a neighbouring cation has the wrong sign (r = +0.11),
and carboxylate proximity gives only r = +0.21 with counter-examples at net +6 sitting at
pKa 3.6. The positions here are empirical, from an exhaustive scan plus a physics solver.

## Selection

Designs passed independent gates before ranking, and the gates **dropped** candidates rather
than down-ranking them — the top-scoring construct of an earlier round was verbatim ubiquitin,
and the top of a 40-design shortlist was 59 % alanine.

- **composition** — model-free; rejects low-complexity and poly-alanine output
- **fold novelty** — FoldSeek + TM-align against CATH50
- **interface** — buried surface area and steric clashes, not contact counts
- **cross-reactivity** — fraction of buried interface area *off* the 9 divergent positions

Ranking follows the competition's own order — pH-selectivity, then cross-reactivity, then
affinity — with one precondition: a design that does not bind cannot show pH-dependent release,
so `tnf2_01` leads as the only candidate carrying both a switch and a credible interface
(923 Å² buried, zero cysteines).

Cross-reactivity, not affinity, did the discriminating: it ranged 0.41–0.89 across variants. One
design with the family's best interface confidence put **59 % of that interface on
human/mouse-divergent residues** and was dropped, along with two others below 0.60.

A measured trade-off worth stating: **each added histidine costs interface confidence** — mean
binder-to-target ipTM 0.235 for singles, 0.180 for doubles, 0.158 for triples. The series spans
that trade-off deliberately rather than picking one point on it.

## What is not established

**No measured binding of any kind.** No expression, SEC, thermal stability, SPR or BLI; every
number above is computational. pH selectivity rests on predicted pKa, and no folding model
represents protonation state. Mouse cross-reactivity rests on a predicted mouse structure and
conserved-epitope coverage, not on measured binding to mouse TNF-α.

Interface confidence is **model-dependent and unresolved**, so a range is reported rather than a
number: AlphaFold2 rated its own optimised designs at i_pTM 0.92 where Boltz-2 rated the same
structures 0.21. No design here was scored by a model that also optimised it. Global ipTM is not
used at all — on a homotrimer target it is inflated by +0.405 on average by the target's own
internal interfaces; all values are binder-to-target chain-pair ipTM, symmetrised and averaged
across diffusion samples rather than taken best-of-n.

One earlier design failed your novelty check at 37.5 % sequence identity while two passed at
38.1 %: the failure was **fold**, not sequence — an immunoglobulin-like β-sandwich at TM 0.893.
That rule shaped this round's selection.
