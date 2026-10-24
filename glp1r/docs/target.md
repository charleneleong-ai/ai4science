# The target: human GLP1R, and why this challenge is not a binder problem

Challenge 3 of the Anthropic × Adaptyv competition asks for a protein that **binds human
GLP1R and switches it on**. Affinity is ranked first and activation second, but activation is
a gate rather than a tiebreak: the organisers state that a binder which does not activate the
receptor does not count. Submissions close **Sunday 25 October 23:59 AoE**, which is Monday 26
October 22:59 in Sydney.

Rules that constrain the design space: de novo and zero-shot only (no starting from or
modifying a known binder, and adequate sequence and structural diversity from known proteins —
known binders may still be used to calibrate filters or train models), 10–250 amino acids, and
one of `protein`, `nanobody`, `scfv`, `fab_kappa`, `fab_lambda`. Up to 20 designs for Tracks 2
and 3, 40 for Track 1.

## The two sites, measured rather than cited

Class B1 receptors are described by a two-domain model: the peptide's C-terminal half lies in a
groove on the extracellular domain, and its N-terminal half drops into the transmembrane bundle
to trigger signalling. Rather than take the groove residues from a review,
[`target_sites.py`](../scripts/target_sites.py) measures them off the active-state structure
(PDB **6X18**, 2.1 Å: chain R = GLP1R 29–423, chain P = GLP-1 7–36, plus Gαs/Gβ/Gγ and Nb35).
Regions follow the UniProt P43220 topology — signal peptide 1–23, extracellular domain 24–139,
TM1 140–164.

**52 receptor residues** lie within 4.5 Å of the peptide, and they split cleanly:

| peptide half | receptor residues | where they are |
| --- | --- | --- |
| **groove**, GLP-1 16–36 | 29 | ECD 19 · ECL1 6 · ECL2 2 · TM1 1 · TM2 1 |
| **activation**, GLP-1 7–13 | 26 | TM1 5 · TM2 4 · TM3 5 · TM5 4 · TM7 4 · ECL2 2 · ECL3 2 |

- ECD groove: `V30 S31 L32 T35 V36 W39 R43 D67 E68 Y69 Y88 L89 P90 W91 R121 L123 E128`, plus
  `P137 E138` at the domain edge.
- Activation site: `L141 L144 Y145 Y148 Y152 R190 V194 K197 L201 F230 M233 Q234 V237 Y241 T298
  N300 W306 I309 R310 I313 D372 R380 L384 E387 L388 T391`.

**The activation site contains no ECD residue at all.** Only three receptor residues (141, 201,
298) are touched by both halves of the peptide, and all three sit at the very top of the bundle.
The two sites are disjoint, and their centroids are **31.6 Å** apart.

Two further numbers describe the activation site itself. Every activation contact lies within
14.0 Å of His7 — the farthest being Y145, at the extracellular end of TM1, so the segment's
contacts are a compact cluster rather than a long smear. And His7 sits **19.2 Å below** the
centroid of the eleven extracellular-loop residues that contact the peptide, which form the rim
of the vestibule it passes through. That 19.2 Å is the insertion depth a design has to cover;
it is a buried position inside the helical bundle, not a surface patch.

## What that geometry rules out

A molecule that does both jobs has to span 31.6 Å and terminate 19.2 Å inside a bundle whose
only access is the loop vestibule. GLP-1 solves this by being a single ~30-residue helix: long
enough to bridge the span, thin enough to thread.

This constrains the format, though less absolutely than the geometry alone suggests. A framework
that presents a flat paratope cannot put anything 19 Å down a channel while staying outside, so
an antibody recognising the groove as a surface patch will bind and not activate. But an
*ultra-long CDR3 loop* is a real theoretical route, and it has been tried at scale: a 10¹⁰
library carrying >150,000 GPCR-binding motifs in HCDR3, median loop length 33–44 residues, gave
13 GLP1R binders and **zero agonists** across five rounds of selection. Affinity maturation of
those got two designs to Emax <2% of exendin-4. So the antibody route is not geometrically
impossible — it is empirically exhausted, which is the stronger argument.

A single-chain `protein` remains the best-supported format, and the inserting element need not
be a terminus at all (see [`strategy.md`](strategy.md) — a free α-amine is not required, and the
de novo class A agonists that worked inserted a *loop*).

This also sharpens the organisers' own stated open question — whether a binder that only fills
the groove switches the receptor on. The groove and the signalling machinery share no residues,
so there is no mechanical path from groove occupancy to bundle rearrangement unless a design
supplies one.

## The affinity/activation split is also an assay split

The two objectives are not measured on the same material:

- **Affinity** is measured against Sino Biological **13944-H18H-DD** — GLP1R 24–432 with an
  N-terminal FLAG and C-terminal polyhistidine tag, purified in DDM/CHS detergent, carrying
  proprietary thermostabilising mutations the supplier does not list.
- **Activation** is measured in whole cells: human cells expressing GLP1R with a secreted
  alkaline phosphatase reporter under a cAMP response element, read as luminescence against
  untreated cells and a reference agonist.

A purified, thermostabilised, detergent-solubilised receptor and a receptor in a live membrane
are not interchangeable targets, and the undisclosed mutations are a real unknown for anything
designed to enter the bundle. This is tracked as an open risk, not a settled fact.

## Not established

- Whether any de novo design has ever agonised a class B GPCR, and whether published designed
  agonists achieved it de novo or by retaining a natural N-terminal motif.
- Whether current structure predictors can model a peptide inserted into a class B TMD pocket
  at all, and whether ipTM means anything for a 7TM target with no membrane.
- Whether the thermostabilising mutations in the assay construct sit in or near the orthosteric
  pocket, and whether they bias it towards the inactive state.
- The ECD-only construct boundaries and whether a soluble-ECD design target is legitimate here.

Each of these is being resolved before any design is generated; none of the numbers above
depends on them.
