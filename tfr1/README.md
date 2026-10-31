# tfr1 — a conditional TfR1 binder, with the epitope chosen by measurement

Adaptyv × Anthropic Protein Design Competition, challenge 04. Three objectives on one binder:
bind human TfR1, cross-react with mouse TfR1, and bind at pH 7.4 while releasing at pH 6.0.
Submissions close **Sunday 1 November 23:59 AoE** (Monday 2 November 22:59 Sydney).

**Status (2026-10-31): target characterised, generation running.** Nothing is measured and no
designs are submitted yet.

## The headline: the recommended epitope serves the lowest-ranked objective

The organisers rank pH-selectivity first, mouse cross-reactivity second, affinity third — and
recommend the apical domain, then note in the next paragraph that mouse differs there at six of
nine residues. Measured here rather than taken on trust (see [`docs/target.md`](docs/target.md)):

| epitope | conserved human/mouse | carboxylates facing it | transferrin clearance |
| --- | --- | --- | --- |
| apical — **the recommendation** | **47%** (VP1u core 22%) | E350, E369 | 46.3 Å |
| 529 — conserved, clear | **90%** (only A529→I) | D525, E533 | 12.3 Å |
| 605 — carboxylate-rich | 83% | **4**: D604 E606 E664 D667 | 6.3 Å |
| 647 — helical | **100%** | D757, E759 | **0.0 Å — it *is* the transferrin site** |

**The most conserved surface on this receptor is the transferrin site.** Conservation and
function coincide, which is the real reason cross-reactive TfR1 binders are hard — the conserved
face is already spoken for. Three epitopes are carried rather than betting the entry on one
reading of the trade-off, and the funnel reports per-epitope survival so the choice is answerable
from data afterwards.

The scan's two binding gates (patch size ≥ 10, transferrin clearance ≥ 6 Å) were **loosened after
seeing the first result**, which is a post-hoc change and is declared as one in
[`docs/target.md`](docs/target.md) with the full sensitivity table — the size gate was landing on
the mode of the patch-size distribution and deciding the campaign by two residues.
`prepare_target.py` now refuses to cut a target whose centre the scan does not mark usable, after
an earlier revision hardcoded all three centres and shipped two that the gate had rejected.

## Calibrated against the organisers' own numbers before use

The scan reproduces two figures stated in the brief, from structure and sequence independently:

- the "nine residues" of the apical patch are exactly the VP1u contacts in **11RN** — 202, 204,
  208, 209, 210, 211, 212, 348, 371 — of which **two are conserved**;
- human/mouse identity over the ectodomain 121–760 comes out at **78.0%**.

Four of the brief's six named divergent residues land exactly; 210 and 211 differ only in where
the gap falls inside the 210–214 indel, which is alignment-ambiguous and leaves both
non-conserved either way.

## The pH switch is designed at the interface this time

The TNF-α campaign ranked its entry on six engineered histidines with pKa in the target window,
then measured the bound-vs-free linkage **Δ = pKa(complex) − pKa(free)** and found its two
best-determined positions **uncoupled on all 25 samples** — they titrate cleanly, and that does
nothing to binding, because an equilibrium identical in both states cancels. All six sat
**9–21 Å from the target, none in contact**.

So the constraint moves upstream: each epitope here was chosen partly for the receptor
carboxylates facing it, which are what a designed histidine can titrate against *at the
interface*. The coupling verifier decides, not the pKa window alone.

## Pipeline

| stage | what runs |
| --- | --- |
| target | [`target_sites.py`](scripts/target_sites.py) measures epitopes, conservation, exposure, glycan and transferrin clearance; [`prepare_target.py`](scripts/prepare_target.py) cuts each target and converts hotspots to the generator's 1-based chain positions, verified per residue |
| generate | BoltzGen, `protein-anything`, 50 designs × 3 epitopes × 2 length bands ([`specs/`](specs)) |
| triage | composition degeneracy, free cysteine, Glu-fraction expression predictor |
| switch | histidine scan at interface positions, then two pKa methods and the coupling check |
| confirm | Chai-1, which shares no weights with the generator |
| novelty | MMseqs2 and Foldseek on the design chain only |

## What is not established

Nothing is measured. The conserved patch at 529 has **no precedent as an epitope** — no structure
shows anything bound there, so unlike the apical site there is no evidence it is druggable, only
that it is conserved, exposed, and clear of both glycans and the transferrin site. The cut targets
expose faces that are interior in the intact ectodomain. Binding is measured on the homodimer, so
avidity may flatter a weak monovalent design. Mouse cross-reactivity is predicted from a sequence
alignment, not a mouse structure.
