# Metal selectivity from occupancy, not from physics (2026-07-14)

The MLIP metal-swap tier is inert: no backbone we have can rank the divalent 3d series
([writeup](2026-07-13-mlip-cannot-rank-metals.md)). That looked like the end of metal discrimination
in touchstone — geometry provably can't do it (a CN4 site is in-range for Ni/Cu/Co alike), so the
`--selectivity` reward had nothing behind it.

It isn't the end. The discrimination was already sitting in the precedent table, unused.

## The signal

For a donor set (say N₂S₂ — His₂/Cys/Met, the type-1 "blue copper" set), ask each metal: **of all the
sites you occupy in the PDB, what fraction use this donor set?** That is
[`motif_enrichment`](../../src/touchstone/geometry/precedent.py), and it recovers coordination
chemistry from pure counting:

| donor set | most characteristic of | the chemistry |
|---|---|---|
| **N₂S₂** (His₂/Cys/Met) | **Cu²⁺** (5.5%) | soft donors — type-1 blue copper |
| **O₆** (all-oxygen) | **Mn²⁺** (6.3%), Cu²⁺ *last* (0.3%) | hard ion prefers hard donors |
| **S₄** (Cys₄) | **Fe** (52%) | iron–sulfur clusters |
| N₃ (His₃) | Cu²⁺ (14.4%), Ni²⁺ (10.8%) | N-rich suits both |

**That is HSAB, with no physics computed.** Hard/soft matching is what MACE and OrbMol failed to
reproduce from first principles; the PDB has it recorded. Costs one table lookup — no GPU.

## Normalising is the whole trick

Raw hit counts are worthless as a metal comparison, because PDB abundance swamps chemistry. For the
type-1 copper donor set, **raw counts say Zn (151 hits) beats Cu (17)** — purely because Zn has ~10×
more sites in the PDB. Dividing by each metal's total presence flips it to Cu, which is the right
answer. Pinned as `test_normalising_by_metal_abundance_is_what_makes_it_right`.

## The smoothing bug that made it lie (caught in review)

The first version smoothed with add-one. That is wrong here, and wrong in a way that produces
confident nonsense: **each metal's denominator is its own total, and those span 60×** (Pt: 53 sites,
Zn: 3105). A constant pseudo-count is therefore not a constant prior — it is a **metal-dependent
floor**:

| metal | total sites | rate implied by **zero** hits |
|---|---|---|
| Pt | 53 | **1.85%** |
| Cu | 309 | 0.32% |
| Zn | 3105 | 0.03% |

A rare metal that has **never been observed** on a donor set outranked an abundant metal with real
precedents. **25 of 85 donor sets named a zero-hit metal as their owner.**

Two fixes: shrink toward the donor set's **pooled background rate** rather than a constant (a thin
metal regresses to *unremarkable*, which is the honest answer when you have no data), and require
`MIN_MOTIF_HITS = 3` actual observations before a metal may be called a donor set's owner. Zero-hit
winners: **25 → 0**.

The lesson is the session's lesson again, in statistical clothing: the failure wasn't a crash, it was
a **plausible number**. An enrichment computed from zero observations is not evidence, and smoothing
made it look like one.

## What it says about the Cu design plan

Something useful, and slightly unwelcome. With Zn²⁺ in the panel, the type-1 copper donor set is only
**1.12×** more Cu-characteristic than Zn-characteristic (Cu 5.5%, Zn 4.9%) — so it is *preferred* but
**not decisive**, and the tier returns `weak`. Against Ni/Co, Cu wins the same donor set by **3.8×**
and trusts.

That is real design guidance: **His₂/Cys/Met buys you Cu-over-Ni, not Cu-over-Zn.** Zinc fingers use
Cys/His too, and the occupancy data knows it. A Cu design that must reject Zn needs something the
type-1 donor set alone does not provide.

## Scope — an occupancy prior, not a binding free energy

This says a donor set is *characteristic* of a metal **in biology**, which is confounded by cellular
abundance and by what evolution happened to need. It **cannot** tell you Cu²⁺ out-competes Ni²⁺ in a
mixed leachate — that is a thermodynamic question, and answering it still needs ligand-field physics
no MLIP we have carries.

So the honest framing of the stack's metal discrimination is:

- **"Does this look like a real Cu site?"** — answered, cheaply, by occupancy. ✅
- **"Will Cu win this site in a mixed feed?"** — still open, still needs DFT. ❌

The first is what an RLVR reward can actually steer on today.
