# Conditional TNF-α binder — method

Adaptyv × Anthropic challenge 02. Three objectives: bind human TNF-α, cross-react with mouse,
and bind at pH 7.4 while releasing at pH 6.0. Ranked pH-selectivity first, then cross-reactivity,
then affinity.

## Epitope choice decides the problem

The three objectives conflict, and the conflict is concentrated in one residue.

The TNFR2 epitope was computed rather than taken from literature: TNF residues within 4.5 Å of
TNFR2 in **3ALQ**, intersected with a human/mouse alignment of the soluble domain (identity
124/157 = 79.0%, reproducing the stated figure). That gives **33 contact positions, 24 conserved
and 9 divergent** (P20H, R31Q, S71D, T72–, H73Y, V85I, T89E, I97V, A111L).

**R31→Q is the trap.** It carries the largest SASA in the entire epitope (193 Å²) and sits beside
R32 in the 29–33 loop that adalimumab and infliximab engage. That is why neither drug is mouse
cross-reactive. Any design gripping that loop wins objective 1 and loses objective 2.

Clustering the conserved, exposed contacts into compact patches gives a better target:

| epitope centre | conserved / divergent SASA | conserved | pH anchor |
| --- | --- | --- | --- |
| **A/91** — A86-92 + B33/B147 | 685 / 72 Å² | **90%** | **K90** |
| B/147 — B20-33, B143-149 | 573 / 68 Å² | 89% | R32 |
| A/115 — 63,67,113,115,143-149 | 378 / 0 Å² | 100% | none |

The **A/91 lobe** was chosen: it is the genuine functional receptor site (preferred per the FAQ),
90% species-conserved by construction, and carries conserved **K90** as a cationic anchor for the
pH switch. One epitope therefore serves all three objectives. Chain C sits 18.7 Å away, so the
target was truncated to the A+B protomer pair (304 residues instead of 471).

## Why the pH mechanism is not "His next to a lysine"

The intuitive design — place a histidine beside K90 so it protonates at pH 6.0 and repels — is
self-defeating. A neighbouring cation **lowers** His pKa, often to ~5.0–5.5, so the histidine
stays neutral at pH 6.0 and the switch never fires. Pairing His with Asp/Glu raises pKa but makes
protonation *stabilising*, which is the wrong direction entirely.

Two mechanisms do work, and both are engineered explicitly:

- **buried** — a His buried *by the interface* pays a desolvation penalty on protonation. This is
  the dominant mechanism in recycling antibodies (Igawa et al. 2010).
- **anchor** — a His held **5–9 Å** from a conserved cation (R32/K90/K128) is repelled when
  protonated. Contacts closer than 5 Å are **rejected**, not rewarded, for the pKa reason above.

Histidine clusters of 2–3 are also emitted, since they raise apparent pKa cooperatively.

**Limitation, stated plainly:** no folding model represents protonation. The `ph_switch` score is
a structural proxy that ranks candidates; it does not predict the pH 6.0 / 7.4 assay.

## Generation

Two tracks, both BoltzGen 0.1.4 conditioned on the A/91 hotspots via `binding_types`.

- **de novo** (`protein-anything`, 60–95 aa). Yield is low: of 8 diagnostic designs, 2 were clean,
  2 collapsed onto **verbatim ubiquitin**, and 3 were poly-Ala backbones the inverse-folder could
  not sequence. Forcing the target into single-sequence mode (`msa: empty`) made collapse markedly
  worse and was removed; it is absent from every canonical BoltzGen example.
- **nanobody** (`nanobody-anything`) grafts CDRs onto 7 real VHH frameworks, holding the framework
  fixed (`structure_groups` visibility 2) and designing only the loops (visibility 0) with
  `design_insertions` for variable CDR length. This **cannot** collapse onto a natural small
  protein, which is why it is the primary track.

## Scoring

All interface quantities derive from per-residue **buried surface area** (unbound SASA minus bound
SASA), never from atom- or residue-pair counts. Pair counting was tried first and was actively
misleading: it scales with interpenetration, so clashing poly-alanine junk out-scored real designs.

| metric | meaning |
| --- | --- |
| `interface_area` | target-side buried area, Å². A real minibinder interface is ~600–1000. |
| `epitope_focus` | share of buried area on the intended conserved lobe |
| `cross_reactivity` | 1 − share of buried area on human/mouse-divergent positions |
| `ph_switch` | interface His weighted by own burial and standoff from a conserved cation |

Designs are hard-gated, not merely down-ranked, on steric clash, low sequence complexity,
interface below 350 Å², and landing off the intended epitope.

## Novelty gating

The competition hard-filters on de novo design, and BoltzGen demonstrably returns memorised
natural proteins, so novelty is a gate rather than a ranking signal. An exact **20-mer** match
against PDB seqres (1.17 M chains) and Swiss-Prot (576 k) is decisive: genuine de novo sequences
essentially never share a 20-mer with a natural protein. It is a memorisation screen, not a
homology search.

The gate is track-aware. Grafted nanobodies deliberately keep a natural framework, so a
whole-sequence screen rejects all of them; CDR1/CDR2 are additionally germline-constrained (the
`YYADSVKG` junction matches nature by construction). The nanobody gate therefore judges **CDR3
only**, at k=12, which is where the specificity and the de novo content sit.

Both gates were calibrated on ground truth before use:

| control | expected | gate |
| --- | --- | --- |
| human ubiquitin | reject | reject (7 hits) |
| human TNF-α | reject | reject (7 hits) |
| synthetic scramble | pass | pass |
| caplacizumab CDR3 | reject | reject |
| designed CDR3 | pass | pass |

## Validation of the engineered mutants

His variants are written back as design files and re-scored through **BoltzGen's own folding
step**, so mutants are judged on the same `iptm` / `interaction_pae` metrics as their parents
rather than by a separate predictor. Every mutant is cross-checked to carry His at exactly the
requested positions, with the structure's sequence matching the submitted sequence — a mismatch
there would mean submitting something other than what was validated.
