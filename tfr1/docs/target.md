# TfR1: the target, and why the recommended epitope is the wrong one

Written before any design was scored, so the epitope choice is a prediction rather than a
rationalisation. Reproduce with [`target_sites.py`](../scripts/target_sites.py) and
[`prepare_target.py`](../scripts/prepare_target.py).

## The objectives are ranked, and the brief's hint serves the lowest one

| rank | objective | measured as |
| --- | --- | --- |
| **1** | pH-selective binding | KD on human TfR1 at pH 7.4, **no detectable binding** at pH 6.0 |
| **2** | mouse cross-reactivity | mouse-to-human KD ratio closest to 1, among designs under 1 µM |
| 3 | affinity | lowest KD on human TfR1 at pH 7.4; binding a functional epitope breaks ties |

The brief recommends the apical domain and, two paragraphs later, says mouse differs there at six
of nine residues — "which is why TfR1 binders are so often species-specific". So the recommended
epitope optimises objective 3 and fights objective 2, which outranks it.

## The apical epitope, measured

The nine residues the brief refers to are exactly the parvovirus B19 VP1u contacts in **11RN**,
which this pipeline reproduces independently: receptor chain A against VP1u chain X at 4.5 Å
gives **202, 204, 208, 209, 210, 211, 212, 348, 371** — nine residues, of which **two are
conserved in mouse (22%)**. Each receptor subunit is engaged by one VP1u copy, so the dimer
presents two independent epitopes.

Aligning P02786 against mouse Q62351 reproduces the brief's overall figure exactly — **78.0%
identity over the ectodomain 121–760** — and four of its six named divergent residues land
exactly. Positions 210 and 211 differ only in where the gap falls inside the 210–214 indel, which
is alignment-ambiguous; both readings make them non-conserved, so nothing downstream changes.

The designed binder 3DS18 (**6WRV**) engages a broader 34-residue footprint that is **62%
conserved** — more conserved than VP1u's, because it spills off the divergent tip. 6WRV carries
the Ser142 allele where our construct has Gly142; 11QC and 11RN match UniProt at every modelled
residue, so both are used in preference and `prepare_target.py` raises on any mismatch.

## Scanning the whole surface instead of accepting the hint

210 of 639 residues are exposed on the **dimer** (RSA ≥ 0.25 — computed on the dimer, since the
dimer interface is not available to bind). Every exposed residue seeds a 12 Å patch, scored for
human/mouse conservation, carboxylate content, and clearance from two things a structure-only
pipeline would otherwise miss:

- the **transferrin site** (619–670), by spatial distance, not sequence overlap;
- the **N-glycans** at N251, N317 and N727. The assay protein is glycosylated and the models are
  not, so a patch under a glycan looks bindable here and is not.

**The most conserved surface on this receptor is the transferrin site.** Patches centred on 644,
647 and 650 are **100% conserved** — and sit 0.0 Å from the transferrin footprint, because that
is what they are. Conservation and function coincide, which is the real reason cross-reactive
TfR1 binders are hard: the conserved surface is the one already spoken for.

### The thresholds were loosened after seeing the scan — declared, not hidden

The first pass gated at patch size ≥ 12 and transferrin clearance ≥ 8 Å and returned **7** usable
patches, the best only 58% conserved. Both gates were then loosened to **size ≥ 10** and
**clearance ≥ 6 Å**, which returns **50**. This is a **post-hoc change**: it was made after
looking at the result, and it is what admits two of the three epitopes below. It is not a
pre-registered filter and must not be read as one.

The reason to loosen rather than accept the first answer is that the size gate was knife-edge —
a 12 Å patch of a 0.25-RSA surface yields ~10–15 exposed members, so the cut landed on the mode
of the distribution and was deciding the campaign by an arbitrary 2 residues. The sensitivity is
therefore reported rather than asserted:

| | clearance 0 Å | 6 Å | 8 Å |
| --- | --- | --- | --- |
| size ≥ 8 | 143 | 113 | 107 |
| **size ≥ 10** | 67 | **50** | 47 |
| size ≥ 12 | 12 | 8 | 7 |
| size ≥ 14 | 1 | 1 | 1 |

Of the 210 scanned patches, rejections break down as 141 on size, 33 on transferrin clearance and
17 on glycan proximity.

### The three carried forward

| epitope | centre | domain | conserved | carboxylates | Tf clearance | glycan |
| --- | --- | --- | --- | --- | --- | --- |
| **conserved** | 529 | protease-like | **90%** — only A529→I diverges | D525, E533 | 12.3 Å | 33.5 Å |
| **carboxylate** | 605 | protease-like | 83% | **4** — D604, E606, E664, D667 | 6.3 Å | 41.3 Å |
| **apical** | 371 | apical | **47%** | E350, E369 | 46.3 Å | 18.6 Å |

529 is the most conserved patch with real clearance from both the transferrin site and the
glycans; 605 is the carboxylate-richest, and subsumes the 604/606 cluster 4–6 Å away, one member
of which (606, 91%) is marginally the most conserved patch in the whole scan; 371 is the brief's
recommendation, carried as the contrast. All three get 50 designs in each of two length bands,
rather than betting the entry on one reading of the trade-off. The funnel reports per-epitope
survival, so the choice is answerable with data afterwards rather than argued.

`prepare_target.py` re-reads this scan at build time and **refuses to cut a target whose centre
is not usable in it**, because the first revision of this campaign hardcoded the three centres and
two of them were in fact excluded by the gate above while being carried forward anyway — a stale
conclusion with no error path. Each hotspot must also be either a member of its patch or an
observed contact; N348 qualifies on the second count, being one of the nine VP1u residues while
sitting just outside a 12 Å patch around 371.

## The pH switch, designed where the last campaign learned it has to be

The TNF-α campaign ranked its entry on six engineered histidines whose pKa sat in the 6.0–7.6
window, then measured the bound-vs-free linkage **Δ = pKa(complex) − pKa(free)** and found the two
best-determined positions **uncoupled on all 25 samples**: they titrate cleanly and that does
nothing to binding, because a protonation equilibrium identical in both states cancels. All six
histidines sat **9–21 Å from the target, none in contact**, which is the most likely reason.

So here the constraint goes in *before* generation rather than being audited after. Each epitope
is selected partly for the receptor carboxylates facing it — D525/E533, the four at 605, E350/E369
— which are what a designed histidine can titrate against at the interface itself. Stage 2 scans
histidine substitutions at interface positions only, and the coupling verifier decides, not the
pKa window alone.

## Not established

Nothing is measured. The conserved patch at 529 has no precedent as an epitope — no structure
shows anything bound there, so unlike the apical site there is no evidence it is druggable, only
that it is conserved, exposed, clear of glycans and clear of the transferrin site. The cut targets
expose faces that are interior in the intact ectodomain; hotspots steer away from them but the
generator is free to ignore that. Binding is measured on the dimer, so avidity may flatter a
weak monovalent design. Mouse cross-reactivity is predicted from a sequence alignment, not a
mouse structure — no experimental mouse TfR1 ectodomain structure is used here.
