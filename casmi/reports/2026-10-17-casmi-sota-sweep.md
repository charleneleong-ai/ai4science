# CASMI 2026 SOTA sweep — closing 0.446 → 0.495 (2026-10-17)

Scope: Kaggle forum (all 101 threads, with comments, pulled through the Kaggle API), the top public notebooks,
and primary-source MS/MS→structure literature from 2024 to 2026. Raw notes:
[`kaggle_discussion.md`](../research_notes/2026-10-17-casmi-sota-sweep/kaggle_discussion.md),
[`methods.md`](../research_notes/2026-10-17-casmi-sota-sweep/methods.md),
[`kaggle_topic_index.txt`](../research_notes/2026-10-17-casmi-sota-sweep/kaggle_topic_index.txt).
All benchmark numbers come from the papers' own authors, not our runs. All forum numbers are single public-LB
submissions or participants' own claims, and none is verified here.

## What changed

- **The public best is now 0.451, not 0.446.** dmitriigluzdov's
  [From Spectra to Structures](https://www.kaggle.com/code/dmitriigluzdov/casmi-26-from-spectra-to-structures)
  (§14) builds on our lineage:
  - v49 0.446 → v51 0.450: depth-7 trees in the second-engine rankers;
  - → v52 **0.451**: standardised [sergiomalv ChEBI/LIPID MAPS](https://www.kaggle.com/datasets/sergiomalv/casmi-chebi-lmsd-structures) set;
  - v53: 8+6 ranker seeds, score not shown.

  So the real gap to the top is 0.044, not 0.049.
- **The gap is close to the noise.** The public LB is 31% of the hidden test (`leaderboardPercentage: 31`, about
  124-132 molecules), so one molecule moving from absent to rank 1 is worth about 0.008. That puts 0.044 at about
  5-6 molecules. A difference between two submissions has a standard error of about 0.016
  ([743254](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/743254)).
  Part of the leaders' lead may not carry over to the private 69%.
- **New host rulings open two routes:**
  - The [ICEBERG PubChem atlas](https://iceberg-ms.mit.edu/) (91.4M PubChem structures, predicted by a
    NIST-trained ICEBERG) "can be used" as model outputs, shipped as a dataset together with the code that queries
    it ([745134](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/745134)).
  - The official [FRIGID](https://github.com/coleygroup/FRIGID) checkpoints are "allowed for prize-eligible
    solutions" ([746430](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/746430)).

  TabPFN, LOTUS, NPAtlas and PubChem PMID/patent/SID counts were also cleared
  ([748882](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/748882)).
- **Treat MassSpecGym numbers with suspicion.** An audit,
  [MassSpecGym in the Wild](https://arxiv.org/abs/2606.19624) (co-authored by Enveda's David Healey), found
  evaluation problems in at least 17 of 26 papers that report MassSpecGym results. The causes are a
  SMILES-format shortcut, a MIST batching bug and leakage. For example, the corrected de novo top-1 for MIST+MolForge
  is 10.7%, against the 31.8% its authors claimed. The audit also reports:
  - PubChem's default ordering alone reaches 50.0% hit@1 on the formula track (author-reported). This matches
    what the popularity prior does in our pipeline.
  - A known formula roughly doubles retrieval: MIST reaches 26.6% hit@1 with a predicted formula and 51.5% with
    the true one.
- **The best forward simulator moved.** [GLACIER v2](https://arxiv.org/abs/2606.29161) (revised 2026-10-05)
  reports, per its authors:
  - MassSpecGym top-1 of 70.0% with contrastive finetuning and 69.7% without;
  - NIST'20 top-1 of 52.5% with finetuning and 38.5% without;
  - about 8 times faster than ICEBERG.

  We should check which GLACIER checkpoint our attached dataset carries. I could not tell from the fork's code.
- **Electron mass.** The hidden test's positive-mode m/z values were corrected by one electron mass. The
  Enveda-180 spectra in train were not, so they need −0.000548 Da to match
  ([748882](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/748882),
  [743395](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/743395)).
  This only matters for tight fragment tolerances, and its effect on our pipeline is unmeasured.

## Kaggle discussion signals

**Test composition.**
- The organisers describe the molecules as "real or potential natural molecules and their analogs"
  ([press release](https://finance.yahoo.com/healthcare/articles/enveda-launching-next-alphafold-contest-120700676.html)).
- The hidden test is timsTOF only, 157-1,159 Da.
- Class 3 (novel) is absent from **both** PubChem and COCONUT
  ([742274](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/742274)).
- The host will not give class shares, but hopes "there are enough novel structures to make it worth working on
  them" ([745758](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/745758)).
  All three classes appear in both the public and private splits.

Participants' estimates of the class shares, all from LB probes and disagreeing:

| Source | Class 1 | Class 2 | Class 3 |
|---|---|---|---|
| prvsiyan / haideptry ([741597](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/741597)) | 16% | 45% | 39% |
| Bull ([748524](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/748524)) | 28% | 20% | 52% |
| Antriksh Arya ([748029](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/748029)) | | | "almost half" |

**What top teams use.** None of the top four (CHU1PC 0.495, pikachu 0.495, Parfionov 0.493, MarvinTMB 0.484)
has posted a method. The indirect signals:
- MarvinTMB reached 0.484 on 2026-10-08, before the public E006 lineage existed. He says "You do not need NIST
  or any additional data", and puts a realistic ceiling at about 0.55
  ([746302](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/746302)).
- hengck23 points to:
  - annotating peaks with fragment formulas before predicting structure;
  - training a forward SMILES→MS/MS model ([744343](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/744343));
  - Enveda's MS2Mol / EnvedaDark work on natural products absent from databases ([745029](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/745029)).
- The host's own [de novo tutorial](https://www.kaggle.com/code/inversion/casmi-denovo-tutorial-notebook)
  (358 votes) is a spectrum→SMILES transformer.

**Measured lessons from the public lineage.**
- Re-ranking gains shrink on the LB:
  - Offline gains under about 0.03 are invisible on the LB.
  - LambdaRank with fragmentation features gained +0.125 on a class-2-like split but only +0.03-0.04 on the LB ([749251](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/749251)).
  - About 98% of errors where the answer is in the list are same-formula isomers ([743254](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/743254)).
- Changes that failed on the 0.440 base, one submission each:
  - RRF with JESTR: 0.418;
  - CFM-ID 4 rescoring against the LOTUS in-silico DB: 0.440 at weight 0.25, 0.422 at weight 1.0;
  - GLACIER per collision-energy step: 0.436;
  - MSEntropy centroiding: 0.440;
  - an independent fingerprint engine merged by RRF: +0.001 (749251).
- Generated candidates have not paid off so far:
  - Inserting them "has not yet paid off on the leaderboard" (743254).
  - Deleting the truth from the DB overstates a generator's reach about 3-fold.
  - A popularity prior is worth +0.50 on library-like validation but about 0 on Enveda-180 molecules without public spectra ([748068](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/748068)).
- ICEBERG's offline gain is about 0 on structures not in MassSpecGym, which suggests it reflects training exposure (748524).
- Environment and probing:
  - An unpinned Docker image cost 0.014 ([749908](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/749908)).
  - A failing assert makes Kaggle score the previous `submission.csv` ([749412](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/749412)).

**What this implies.** The public lineage is near the ceiling for re-ranking its own pool. What it leaves on the
table:
- class 3, which no pool can reach;
- class-2 isomers, where structure-aware evidence such as forward models and fragments is the only lever left.

The leaders' 0.04-0.05 is most plausibly a few extra class-2 isomer wins plus a few class-3 hits. That is an
inference from the evidence above, not something anyone has confirmed.

## Methods (table)

H@k means hit@k in %. "Mass" and "formula" are MassSpecGym's two candidate tracks. De novo numbers assume the
true formula is given. More detail and the URLs are in
[`methods.md`](../research_notes/2026-10-17-casmi-sota-sweep/methods.md).

| Method | Type | Author-reported numbers | Weights / licence / eligibility | T4, about 400 molecules | Plug-in point |
|---|---|---|---|---|---|
| [GLACIER](https://arxiv.org/abs/2606.29161) | forward simulator | MSG top-1 70.0% (contrastive), NIST'20 52.5% | ms-pred MIT; MSG checkpoint public; host-approved | 369 spectra/s on A5000, so the full pool is feasible (estimate) | already used for top-40 same-formula reordering → also score the **whole shortlist** as a ranker feature, and check the checkpoint |
| [ICEBERG](https://www.biorxiv.org/content/10.1101/2025.05.28.656653v1) / [PubChem atlas](https://iceberg-ms.mit.edu/) | forward simulator / precomputed library | NIST'20 [M+H]+ top-1 40.0%; scaffold split 33.5% | the atlas comes from a NIST'23-trained model; **host allows its outputs** (745134). Atlas spectra are Orbitrap NCE, PubChem only, [M+H]+/[M−H]− | no inference needed; offline download and binning only | precomputed spectral library for PubChem candidates → ranker feature |
| [FRIGID](https://proceedings.mlr.press/v306/bohde26a.html) | de novo (MIST fingerprint → masked diffusion LM + ICEBERG refinement) | MSG top-1/10 18.3/22.0% (valid per the audit); NPLIB1 25.0/33.4% | code Apache-2.0; Zenodo weights CC BY-NC but **host-approved** (746430); MSG/NPLIB1 data | 6.6 s/spectrum on RTX 6000 Ada, about 3-4.5 h on T4 for all; gate to low-confidence molecules | class-3 candidate generator; needs formulas → take the top formulas from the pool |
| [MS-GPT](https://arxiv.org/abs/2607.23607) | de novo (SAFE-GPT conditioned on fingerprint + formula) | MSG top-1/10 23.9/28.7%; NPLIB1 29.8/41.1% (not audited) | code Apache-2.0; weight licence unstated; bundles FRIGID's MIST encoders (third-party mirror; reproducibility caveat, 746430) | 100 samples per spectrum; about 1-2 h (estimate) | alternative to FRIGID; same plug-in point |
| [FLARE](https://pmc.ncbi.nlm.nih.gov/articles/PMC12873900/) | peak↔atom contrastive | MSG mass 43.2/75.6/92.9, formula 22.7/50.0/75.2 | MIT, trained on MSG | cheap | ranker feature; same family as JESTR, which **lost** on the LB (0.418), so low priority |
| [MIST](https://github.com/samgoldman97/mist) / MIST-CF | spectrum→fingerprint / formula | MIST vFRIGID 53.8 formula H@1 (audit table) | MIT; bundled weights' data unclear; batching bug | 0.6 s/spectrum | second fingerprint engine; the evidence says independent engines add about +0.001 |
| [DreaMS](https://github.com/pluskal-lab/DreaMS) | SSL spectrum encoder | no clean MSG retrieval number; frozen alignment drops to 6.7 H@1 (audit) | MIT, GeMS (no NIST); approved | cheap | backbone for a fine-tuned fingerprint head |
| JESTR / MVP / SpecBridge / GLMR | contrastive | MSG claims, most likely inflated by the format shortcut | mixed | cheap | skip (JESTR measured −0.022 on the LB) |
| DiffMS / MADGEN / MS-BART / MBGen / MSNovelist / FOAM | de novo | MSG top-1 0-7.6% | mostly MIT | DiffMS 131 s/spectrum; FOAM about 3 min/spectrum | skip |
| CSI:FingerID / SIRIUS 6 / COSMIC | fingerprint + DB | CASMI 2022: SIRIUS team 26-29% top-1, formula 94-95% | needs a web-service account | not offline | skip (only fragmentation trees run locally) |

## Ranked interventions

Ranking is expected gain × feasibility within about a week, at 5 submissions/day and about 24 GPU-h. The
expected gains are rough estimates from the evidence cited, not measurements.

| # | Intervention | Evidence | Expected public-LB gain | Cheapest first test |
|---|---|---|---|---|
| 0 | **Catch up to 0.451**: depth-7 second-engine trees, the sergiomalv ChEBI/LMSD set, 8+6 seeds, pinned image | dmitriigluzdov §14 (v51 0.450, v52 0.451) | +0.005 (already measured publicly) | 1 submission of the fork with the three changes; check the image pin first |
| 1 | **FRIGID as a gated class-3 generator.** Run only where retrieval confidence is low (max library similarity < ~0.7 and a small FPNet margin). Condition on the top 2-3 formulas of the fused list. Insert its 1-2 best candidates at **ranks 3-5**, never rank 1, deduplicated by the metric key | Class 3 is 39-52% of the test (probes) and no pool reaches it. FRIGID is the strongest de novo result that holds up under the audit (NPLIB1 25% top-1, natural products); host-approved. Counter-evidence: generated inserts have not paid off for others (743254) | Back-of-envelope: 45% class 3 × 10-20% top-5 hit rate × about 0.3 reciprocal rank ≈ **+0.013-0.027**. The cost is a one-rank shift for gated molecules only | Offline on Kaggle T4: 100 held-out Enveda-180 structures with no public spectra (Aryan's S2 recipe), formulas taken from our pipeline. Measure top-5 exact match with **predicted** formulas and runtime. Submit only if top-5 is ≥10% |
| 2 | **GLACIER over the whole shortlist as a ranker feature**, not just top-40 same-formula reordering. Use the contrastive-finetuned checkpoint, score every collision energy and take the max, cover [M−H]− with ICEBERG | 98% of the remaining errors are same-formula isomers (743254). GLACIER is the strongest retrieval model reported (70% MSG top-1). Earlier GLACIER steps moved the LB (v18 +0.007, v40 +0.007, E006 +0.018, though that step also added metric-identity fusion). Caveat: the gain on structures outside MassSpecGym may be near zero (748524) | +0.005-0.015 | Zero submissions: compare the attached checkpoint with the ms-pred release. Then re-score the isomer panel (Enveda-180 molecules without public spectra, **not in MassSpecGym**) with top-40 versus a 200-candidate shortlist |
| 3 | **ICEBERG atlas as a NIST-grade in-silico library** for PubChem candidates (class 2): download atlas spectra by formula for the formulas the pool can reach, bin them, and ship the dataset plus the query script | Host-approved (745134). The NIST'20-trained ICEBERG reports top-1 40% versus weaker MSG-trained models. Public notebooks do not use it, so it is a plausible lever for the leaders (**unverified**). Caveats: Orbitrap NCE versus timsTOF eV; PubChem only (COCONUT-only structures missing); data engineering | +0.005-0.015 for class 2 (speculative) | Pull atlas `.mgf` files for the formulas of 100 held-out panel molecules and compare same-formula ranks with MSG-GLACIER. The offline cost is download time only |
| 4 | **Reach check on new reference spectra (class 1).** Take newer GNPS/MassBank/MoNA releases plus the 816 Enveda-180 Zenodo structures missing from train (host-allowed, 745148), and count InChIKey14s that are new versus train | A class-1 hit scores about 0.93 MRR versus about 0.3 otherwise, so each newly covered test molecule is worth about +0.005 on the LB | Unknown; probably small, because train compiles the same libraries | Zero submissions: dedupe the releases against train and count new structures in 157-1,159 Da. Submit only if the count is large |
| 5 | **Fragment-formula / subformula annotation features** (hengck23, Adarsh). Annotate peaks with subformulas of the candidate formula, then score candidates by the explained intensity of 1-2-cut fragments relative to the best candidate in the molecule. Add the −0.000548 Da positive-mode fix for Enveda-180 references | Aryan: the relative fragment-gap feature was his strongest (S2 0.502 → 0.627), with an LB effect of about +0.03-0.04 across his steps. The NeckBeard tweaks held out at +0.023 panel. Our pipeline already has MetFrag-lite, so the extra gain is partly taken | +0.003-0.008 | Add a "fragment gap to best candidate" feature and the tolerance tweaks to the engine-2 ranker, retrain with 8 seeds, compare on the identity-disjoint folds; 1 submission |

The order of work follows from this: #0 now; the zero-submission checks for #2 and #4 in parallel (CPU and
download time); the FRIGID offline probe (#1) on T4 time; then one submission per verified change.

## What would change the plan

- **FRIGID's top-5 with predicted formulas is under about 5%** on the Enveda-180 no-public-spectra panel. Then
  class 3 is out of reach in a week: drop #1 and put effort into #2/#3.
- **The attached GLACIER is already the contrastive v2 checkpoint, and shortlist expansion moves nothing on the
  non-MassSpecGym panel.** Then drop #2.
- **A top team posts a method, or the host releases class-stratified hints.** Then refresh this sweep.
- **The gap is mostly public-LB noise.** If our own CV puts 0.446 → 0.451 changes inside seed noise, the
  priority shifts. Robust final-submission selection matters more than chasing 0.495: the private LB is 69% of
  the test, and the leaders' 0.04 equals about 2.5 standard errors.
- **Eligibility.** The atlas route needs the query code and the predicted spectra shipped as a dataset. If the
  host tightens the rules on NIST-trained outputs, #3 becomes medal-only.

Not established:
- the class shares;
- whether any leader uses generation or the atlas;
- FRIGID and MS-GPT performance on timsTOF CID spectra with predicted formulas;
- all T4 runtimes, which are estimates;
- MS-GPT's weight licence.
