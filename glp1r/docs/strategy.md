# Strategy: what to build, and what the literature says will happen

Written 2026-10-17, before any design was generated, so that the filters and the portfolio split
are on record ahead of the results. Target geometry is in [`target.md`](target.md).

## 1. The obvious entry is a published negative control

Six independent campaigns have put a protein in the GLP1R extracellular-domain groove. **None
produced a full agonist.**

| ligand | epitope | function |
| --- | --- | --- |
| exendin-4(9-39) | ECD groove, IC50 6 nM at the isolated ECD | clinical **antagonist** |
| Fab 3F52 (5E94) | ECD groove, ~958 Å² buried | **neutral antagonist**, no agonism to 1 µM |
| Fab 7F38 (6LN2) | ECD, off-groove | **silent** — neither blocks nor potentiates |
| Glp1R0017 (phage scFv→IgG1) | ECD, IC50 5.2 nM | **antagonist** |
| dGl1_024 / mGl1_008 (de novo) | ECD, hotspot L95 | **antagonists**, 61 and 39 nM |
| TB-222-040 / -089 (HCDR3 helix mimics) | ECD (inferred) | **partial, Emax <2% of exendin-4** |

The strongest single result is a 10¹⁰-diversity antibody library carrying >150,000 GPCR-binding
motifs grafted into HCDR3: five rounds of selection on GLP1R, 13 specific binders, **zero
agonists**, stated verbatim by the authors. In that same study agonism appeared only when native
GLP-1 was physically fused to the light-chain N-terminus — which is how every marketed
antibody–GLP-1 agonist works, and is not de novo.

And the de novo version has been run. Muratspahić et al., *Nature* 656:1044–1053 (2026) took
**10,080 designs → 49 expressed → 2 antagonists** against the soluble GLP1R ECD, read out in a
CRE-luciferase line — the same assay family as this challenge's CRE-SEAP. Their stated rationale
was that ECD binding should sterically block the peptide and *produce* antagonism. No class B
agonist appears anywhere in that paper; their only agonists were class A (MRGPRX1), and the
reason they give is that agonism needs "conformational selectivity between active and inactive
orthosteric states".

**So a high ECD-groove affinity score is not progress toward the activation gate.** It is the
single best-evidenced claim available about this target.

## 2. The switch is in the bundle, and the pharmacophore is rankable

An 11-mer corresponding to GLP-1 7–17, carrying **no ECD-binding helix at all**, is a 0.5 nM
full agonist (optimised 0.01 nM). The groove is the affinity trap; the bundle is the switch.

Sequence identity to the natural motif is **neither necessary** — a novel 9-residue extension,
`ELVDNAVGG`, gave a 0.85 nM full agonist at Emax 102%, and saturation mutagenesis of GLP-1's
first five residues (~1700 variants) returned 27 functional agonists — **nor sufficient**: two
designs retaining the complete His7–Phe12 motif were both completely inactive, and the 10¹⁰
antibody library that grafted the GLP-1 helix into HCDR3 produced zero agonists from 13 binders.

So the gate must be geometric and chemical, never sequence identity. From the mutational and
structural record, in priority order — this is the design brief:

| rank | feature | receptor partner | cost of losing it |
| --- | --- | --- | --- |
| 1 | **anionic group three residues in** from the inserted end | bidentate salt bridge to **R190** (2.78/3.14 Å) + H-bond **Y152** (2.76 Å) | the one feature >80% of randomised-N-terminus agonists converged on |
| 2 | **aromatic** ring | **Y148** T-stack, in the L141/L144/L388 slot (Phe12 site) | Ala here is >1000-fold in the 11-mer scan |
| 3 | **glycine-sized** atom | **W306** (Gly10 slot) | Gly10→Ala abolishes activity |
| 4 | **imidazole-like** H-bond donor | **Q234**, inside the V237/W306/I309/I313 cage (His7 site) | 10–200-fold; not absolute |
| 5 | small hydrophobe | **E387/L388** (Ala8 methyl) | nearly free |
| 6 | ECL3 contacts **D372/E373/R380** | closes TM6–ECL3–TM7 | separates *full* from *partial* efficacy and unbiased from biased |

**A free α-amine is not required, and the inserting element need not be a terminus.** N-acetyl
GLP-1 costs 89-fold in affinity but only 3.5-fold in cAMP potency; the 10-mer agonist in 5NX2
has no α-amino group at all, its first residue being a histamine/gem-dimethyl-malonyl surrogate
whose *carboxyl* forms the peptide bond. At 2.1 Å the His7 α-amine makes no direct hydrogen bond
to the receptor — its nearest polar partner is 5.3 Å away — so the widely repeated claim that it
salt-bridges E387 or E364 traces to models, not coordinates. And the de novo class A agonists
that worked inserted a **loop**, not a terminus: maxadilan agonises PAC1R with its central turn
while its helix ends hold the ECD.

This widens the design space considerably. A protruding loop on a folded mini-protein is far
easier to build than a solvated free terminus threaded down a channel.

**The template to copy is NK1R.** Free RFdiffusion against the **active state** gave 20 agonists
from 63 expressed designs (EC50 1–231 nM); the other generator, which produced no
substance-P-like motif, gave **0 active from 20**. What separated them was independent
convergence on the natural ligand's activating loop geometry — without its sequence. That is
exactly the filter described in §4, and it is the only thing in the published record that
actually discriminated agonists.

**Hard constraint:** at GLP1R the ECD anchor is *not* optional. The isolated TMD shows no
activation up to 10 µM, even with hormone tethered. Both halves are load-bearing, which is why
the 31.6 Å span is a requirement rather than an inconvenience.

## 3. How to get the geometry without copying the ligand

The rules require de novo and zero-shot: no starting from a known binder, no modifying one,
adequate sequence and structural diversity from known proteins. They explicitly permit known
binders for **calibrating filters or training models**.

The published recipe for a deep pocket is motif scaffolding — place a short peptide motif on the
pocket hotspots, freeze it, build the protein around it. Seeding that motif with GLP-1's own
backbone would be starting from the natural ligand, so this campaign does not do it. Instead:

- **Generate** from the *receptor* pocket hotspots only, with no ligand template.
- **Select** on recapitulation of the activating geometry, with the natural ligand used purely
  as the filter's reference.

That keeps the generator zero-shot while still applying the one lesson that separated working
designs from dead ones in the published record, and it puts the natural ligand on the
filter-calibration side of the rule, where it is allowed.

## 4. Scoring: almost everything available is miscalibrated here

This is the part most likely to sink the entry, and it is measurable.

- **Do not gate on ipTM.** Against 601 *wet-lab-tested* designs from Adaptyv's own earlier
  rounds, ROC AUC for binding was ipTM **0.64**, pLDDT 0.66, alntmscore 0.66 — and the best
  metric of all was sequence identity to the PDB at 0.68. The authors found no significant
  correlation between Kd and iPAE, ipTM or ESM2 likelihood, and call iPAE and ipTM
  "insufficiently predictive". The two highest-AUROC metrics being similarity-to-known-protein
  measures is memorisation surfacing in a wet-lab readout.
- **Use inter-chain PAE and binder pLDDT instead.** On 113 GPCR–peptide complexes (including 5
  GLP1R), Spearman ρ against DockQ: inter-chain PAE −0.894 (RF3), −0.829 (AF2 initial guess),
  −0.56 (Boltz-2); peptide pLDDT 0.76 for Boltz-2 (0.81 best-seed).
- **Expect mini-protein scale to be much worse than peptide scale.** On ligands >50 residues —
  our size class — **35.5%** of Boltz-2, 82.1% of RF3 and 100% of AF2-initial-guess predictions
  scored DockQ < 0.23. And post-training-cutoff, high-quality share collapses: Boltz-2
  20.8% → 6.4%. A novel design is post-cutoff by construction.
- **On membrane proteins Boltz-2's own confidence is nearly uninformative** (Pearson r = −0.095
  against ligand RMSD). The 7TM bundle itself predicts well without lipid — TM Cα ~1.0–1.4 Å —
  but what fails is the pocket *contents*, and the failure is in sampling, not ranking.
- **Boltz-2's affinity head cannot score this binder at all.** It requires the binder chain to be
  a small-molecule ligand, caps at 128 atoms and was trained at ≤56; a 60-aa mini-protein is
  ~450 heavy atoms. Zero time goes here.
- **Never score a design with the model that generated it.** Pairing a generator with its own
  family's validator measurably inflates results. BoltzGen generates *and* ranks with Boltz-2
  weights, so an independent co-folding model is required, not optional.

### What that implies

1. **Calibrate before trusting any score.** Co-fold a panel whose functions are known and
   opposite — GLP-1 and exendin-4 (full agonists), exendin-P5 (full agonist, novel N-terminus),
   GLP-1(9-36) (weak partial, EC50 363 nM at Emax 85%), exendin-4(9-39) and Fab 3F52
   (antagonists) — and ask whether any computable score separates agonists from antagonists and
   ranks the partial between them. If none does, say so and do not rank on it.
2. **The cheapest high-value filter is the geometric gate written by hand**: is the design's free
   N-terminus actually buried inside the bundle, at depth, with the R190 salt bridge — and not in
   the ECD groove, and not on the intracellular face? 21.7% of Boltz GPCR predictions land in the
   *wrong pocket* with plausible confidence, and templates have biased peptides onto the
   intracellular side (one case: DockQ 0.02 with templates vs 0.60 without). This must be
   checked on every design, every time.
3. **Cross-model pose consistency beats any single confidence score.** Agreement across
   independent predictors reached 91% high-quality selection where single scores could not.

## 5. The portfolio

The assay is a two-stage funnel — activation is measured only on confirmed binders — so every
design needs ECD-groove affinity to survive stage 1 *and* a bundle-reaching element to pass
stage 2. One architecture, with insertion depth as the deliberately varied axis. The ECD anchor
is common to all tiers, so the arm is additive rather than traded against affinity.

| tier | what it is | why it is in the entry |
| --- | --- | --- |
| **groove-only** | ECD anchor, no arm | affinity anchors, and the negative control that makes the depth series interpretable |
| **shallow shelf** | ECD + upper TM1/TM2/ECL1 | best risk-adjusted shot at the gate: GLP-1(9-36) binds exactly this shallow site and still reaches EC50 363 nM at Emax 85%, and a CRE reporter carries ~3 orders of magnitude of receptor reserve |
| **deep** | free N-terminus into the His7 pocket, R190 salt bridge | the only route with a mechanism for full agonism |

The split is recorded here before results exist. This makes the entry a pre-registered test of
the organisers' own stated open question — whether a groove-filling binder switches the receptor
on — rather than a bet dressed up as a prediction.

## 5b. The first gate is selection, not the wet lab

Terms §5.4 separates two numbers that are easy to conflate. The *submission* cap is 40 for
Track 1 and 20 for Tracks 2–3. The *testing* allocation is different: Track 1 gets ~15 designs
tested per problem per team, while **Tracks 2 and 3 share a single 384-well plate per problem,
at the Sponsors' discretion, with no guaranteed allocation**, from a pool of roughly 375 designs
per track. Selection is made by a Claude-based workflow whose prompt is deliberately unpublished
"to avoid teams over-optimizing metrics", judging predicted design quality, novelty **and method
novelty** — and it reads the submitted information. Attempting to influence it by injection is
disqualifying.

Three consequences:

- **The methods narrative is load-bearing**, not documentation hygiene. A pre-registered depth
  gradient that tests the organisers' own stated open question is method-novel in a field where
  the public competitor repos are uniformly BindCraft on soluble targets.
- **Diversity across the 20 matters.** A pooled selector choosing from ~375 will not take 20
  near-identical helices.
- **The novelty scorer is a late-stage killer.** Published thresholds: "de novo" requires
  sequence identity ≤30% *and* less than moderate structural similarity, scored by MMseqs2
  against SwissProt/PDB/patents/PLAbDab and Foldseek/TM-align. Participants report designs
  rejected at 2/4. Every slot gets checked before it is spent.

## 5c. Two measured facts about the assay construct

**Peptide agonists bind the detergent-solubilised receptor only weakly.** Two independent
vendor datasheets for the same UniProt 24–432 range in DDM/CHS report dulaglutide at **1.91 µM**
and a GLP1R agonist at **2.74 µM**, against **69 nM** for dulaglutide on the challenge's own
construct — and the isolated ECD binds exendin-4 at **22–58 nM**. The ECD groove is therefore
the stronger binding epitope on this material, not the bundle, which pulls affinity and
activation in opposite directions exactly as §1 implies.

**The ECD is wild-type in every published thermostabilised construct.** Across all three
disclosed lineages — the ten-mutation inactive-state set, and both Heptares StaRs — not one
mutation falls in residues 24–145. Several *do* line the orthosteric pocket (T298A contacts
GLP-1 Ser17 at 2.46 Å; M233C, E387D, Q211A, D215R), and the inactive-state set deliberately
locks TM6 via an I317C–G361C disulfide, measured at 2.05 Å. The supplier discloses none of its
own mutations — one boilerplate sentence, no count, no positions, no Tm data, and no native-peptide
binding figure on the datasheet at all. Treat the composition as unknown.

## 5d. ECD anchor hotspots, from three independent sources

- **L95** — the hotspot the two published de novo ECD antagonists were designed against.
- **W33** and the **L32–R43 helix** — Fab 3F52's epitope core. W33 is also the primate-specificity
  residue for both oral agonists. Worth noting W33 sits at 5.20 Å from GLP-1 in 6X18, just
  outside the 4.5 Å cutoff used in [`target.md`](target.md), so it appears as a validated hotspot
  in the literature but not in our own contact set — a case where measurement and literature
  must both be consulted.
- The buried-area ranking within the groove: core hydrophobic **L32, W39, W91, Y69**; rim
  charges **E68, R121, E128, P90, T35**.

## 6. Honest expectations

- **No de novo class B1 agonist exists.** Not GLP1R, PTH1R, GCGR, GIPR, PAC1R, CRF1R or CGRPR.
  De novo GPCR agonism exists only at class A receptors with shallower pockets, and the strongest
  is EC50 390 nM.
- **Deep-pocket binding is rare.** Orthosteric-pocket epitopes gave ~5–10% binder hit rates for a
  closed frontier model where surface epitopes gave 48–100%. The one de novo deep-7TM success
  needed ~26,000 designs for 2 hits, and it was an antagonist with explicitly "no agonism".
  Published deep funnels run 13k–100k designs; that is one to two orders of magnitude beyond
  nine days on one A100.
- **Surface-target hit rates are the realistic anchor.** Adaptyv's own earlier rounds ran 2.5%
  (5/201) and ~13%. The ECD is a clean soluble groove, which is the regime the tools are built for.
- So the realistic win condition is **landing binders plus a well-posed depth experiment**, with
  a genuine but low chance of a weak partial agonist. It is not "first de novo class B agonist".

## 7. Not established

- Whether BoltzGen's pocket conditioning reaches a buried class B site at all — there is no GPCR
  or membrane precedent for it upstream, so this is a genuine experiment and gets a small pilot
  before any budget.
- Whether the assay construct's undisclosed thermostabilising mutations sit in the orthosteric
  pocket. Circumstantially, 5NX2 — a thermostabilised GLP1R — spans UniProt 24–432 exactly, which
  is the range the supplier states, with 12 flagged mutations. Not proof.
- Any class-B1-stratified co-folding accuracy number. The 113-complex benchmark contains 14 class
  B1 receptors and 5 GLP1R complexes but reports only class composition, not stratified DockQ.
- Whether any structural proxy predicts agonism. Every activation-state classifier published is
  trained and tested on *experimental* structures; none is validated on predicted ones, and
  partial agonists were excluded from the best one's training set. No TM6/NPxxY/G6.50 metric has
  been validated as a *predictor* of a novel ligand's efficacy — they are post-hoc descriptions.

## 8. Tool state on the box

Verified by direct inspection and smoke test, 2026-10-16/17.

Ready, no download: **BoltzGen 0.1.4** (the only installed tool that can express a fixed literal
N-terminal segment plus a designed remainder, and crop the target to a sphere around a chosen
residue); **Boltz-2 2.2.1** with weights; **BindCraft2 1.0.3** — AF2 params resolve with
`BINDCRAFT_AF2_PARAMS=~/BindCraft/params`, confirmed, no download; **BindCraft v1**, which is
broken by default because a user-site numpy shadows the env and is fixed by `PYTHONNOUSERSITE=1`,
confirmed; **RFdiffusion2** backbone generation via its 13.6 GB container.

Caveats carried forward: installed BoltzGen is 0.1.4 against upstream 0.3.2, so the spec grammar
needs re-validating before a large run; BoltzGen indexes by mmCIF `label_seq_id`, not author
numbering; and a prior campaign on this stack recorded that omitting `total_len` makes it spin
forever. BindCraft2 ships `termini_accessible`, which points termini *away* from the target — the
opposite of what the deep tier needs — so it cannot express the insertion constraint and is
scoped to the groove tier.

Downloading: RFdiffusion2's missing sequence-design and validation containers (14.6 GB), which
were impossible at 9.5 GB free and fit after reclaiming 48 GB. The validation container carries
**Chai-1**, which also supplies the independent co-folding model that the self-confirmation
problem in §4 requires.
