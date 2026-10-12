# Experiment card — conditional TNF-α binder

Structured after [*A generated protein is not a result*](https://rewirebio.io/blog/a-generated-protein-is-not-a-result/).
A generated sequence is a proposal that passed computational checks; a result is a measured
biological endpoint. Nothing here is a result.

## 1. Target, frozen

| item | value |
| --- | --- |
| target | human TNF-α soluble homotrimer, UniProt **P01375-1**, residues 77–233 (157 aa mature) |
| structure | **1TNF** chains A/B/C, 2.6 Å. Residues 1–5 disordered; modelled 6–157 |
| receptor epitope | TNF residues within 4.5 Å of TNFR2 in **3ALQ** — 33 positions |
| ortholog | mouse TNF-α, UniProt **P06804**; trimer predicted (no experimental structure used) |
| epitope targeted | conserved pocket 113, 115, 143, 144, 146, 147 (and a 86–92 lobe variant) |
| coldspots | the 9 human/mouse divergent positions: 20, 31, 71, 72, 73, 85, 89, 97, 111 |

A **1TNF vs UniProt discrepancy** is recorded rather than smoothed: position 143 is Leu in the
crystal structure and Asp in P01375.

## 2. The generated object

Single-chain de novo miniproteins, 61–110 aa, `molecule_class = single_chain`. Histidines are
then engineered at interface positions to create the pH switch.

## 3. Stack, exactly

| role | tool | version / checkpoint |
| --- | --- | --- |
| design (sampling) | BoltzGen | 0.1.4, `boltzgen1_adherence` + diverse checkpoints |
| design (optimisation) | BindCraft | 1.5.3, AF2 `params_2022-12-06`, `default_4stage_multimer_hardtarget` |
| design (dual-target) | BindCraft2 | repo `bindcraft==1.0.3`, jax 0.11.2 |
| co-folding / validation | Boltz-2 | 2.2.1, `boltz2_conf.ckpt` (affinity head is ligand-only, unused) |
| side-chain building | PyRosetta | 2026.29 release |
| pKa prediction | PROPKA | 3.5.1 |
| sequence novelty | MMseqs2 | vs Swiss-Prot + PDB seqres, snapshot 2026-10-10 |
| fold novelty | FoldSeek + TM-align | CATH50, snapshot 2026-10-11 — **conservative vs Proteinbase**, which uses ESMFold + domain segmentation over PDB/patent/antibody sets |

Terms: BoltzGen/Boltz MIT; BindCraft MIT; **PyRosetta requires a licence — free for
non-commercial academic use, paid otherwise**; AlphaFold params CC BY 4.0.

## 4. Sampling budget and the funnel, with denominators

Failed generations are in the denominator.

| stage | n | note |
| --- | --- | --- |
| BoltzGen designs generated | 100 | across 4 screen configs + 2 earlier runs |
| campaign constructs (parents + His variants) | 57 | |
| passed sequence-novelty (exact 20-mer) | 53 | 4 memorised — 2 ubiquitin, 2 VHH framework |
| passed composition gate | 7 | 46 lost to Ala/Gly degeneracy |
| passed structural quality (clash / interface / epitope) | **3** | |
| submitted | 3 | |
| passed Proteinbase auto-novelty | **2** | 1 failed on fold similarity |

**Yield from generated to novelty-passing: 2/57 ≈ 3.5%.**

Histidine scan: 200 single substitutions → 18 ideal (pKa 6.3–7.2), 17 usable (6.0–7.6).

## 5. Oracles, thresholds, and where they disagree

Interface confidence is **not** reported as a point estimate, because the models disagree
on identical structures:

| structure | AF2 (BindCraft) | Boltz-2 | BoltzGen folding |
| --- | --- | --- | --- |
| BindCraft design | **0.92** | **0.21** | — |
| BoltzGen designs | — | 0.44–0.64 | ~0.20 |

Two artefacts were found and corrected, and both would have inverted conclusions:
- **Global ipTM is unusable on a multi-chain target.** The TNF trimer's own internal interfaces
  score ~0.88 and swamp the binder's. Moving dimer→trimer lifted global ipTM 0.20→0.52 while
  the binder-to-target interface stayed flat. All numbers here are **binder-to-target only**.
- **Best-of-n inflates.** Metrics are stored per diffusion sample; `best` and `mean` are both
  reported and ranking uses **mean**.

## 6. Filters and ranking

Hard gates (artefacts, not weak candidates): memorised sequence; single residue ≥30% or <12
distinct residues or Ala+Gly >45%; steric clash; interface <350 Å²; off-epitope; fold-novelty
score <3.

Ranking: pH switch, then binder-to-target ipTM (mean), then cross-reactivity — following the
competition's own order, and its statement that a weak-but-selective binder outranks a strong
non-selective one.

**Declared deviation:** these filters were **not** all fixed before seeing results. The
composition gate was added after poly-alanine appeared at the top of a shortlist, and the fold
gate after a submitted design failed. Both were responses to observed failures, so the funnel
above is not a clean pre-registered statistic.

## 6b. Second submission, and seven more declared deviations

Submission 2 is a measured histidine series on the two chassis that passed **both** the fold gate
and the composition gate — of 16 folded parents, only `screen_pocket_helix_0` and
`screen_pocket_free_1` passed both. Thirteen entries vary the switch on `helix_0`, two carry
`free_1`, and five of the 20 slots go unused. The concentration on one chassis is a result,
not a preference: the funnel left two.

The seven unused slots were held for the BindCraft2 campaign, which produced nothing in three
rounds (§ declared deviations). Rather than leave them empty on a guess, the eight remaining
three-histidine combinations on `helix_0` that hold all three pKas in window on their design
structure were folded against human and mouse TNF-α and scored on the same gates as the rest
of the set. **One of the eight survived** and is submitted as `tnf2_13`; the full funnel,
including each rejection, is in `docs/slotfill-funnel.csv`. Six failed on structure (clashes,
or an interface that moved off the epitope entirely) and one on mouse retention at 56 %.

A second round of eight closed the one assumption the first round left standing. Under the
clash gate as originally set it returned **zero survivors** (`docs/a38-funnel.csv`); under the
gate recalibrated to the submitted set's own tolerance it returns **one**, `tnf2_14`, which has
the largest guaranteed swing of any design in the submission. Six A38H triples and two R21H doubles had been
dropped earlier without being folded, on the grounds that `helix0_E38H` as a *single* put 59 %
of its interface on human/mouse-divergent residues — an inference from the conserved-epitope
proxy this campaign later disproved. The exclusion was correct and **the stated reason was
wrong**: cross-reactivity was never the problem, with five of the eight at or above 0.70 and
the best at 0.87. Seven fail on structure, and the one that is structurally clean
(`helix0_R21H_E38H_E67H`, 0 clashes, epitope focus 0.48) fails on robustness at 1 of 3
histidines. Across both rounds, **3 of 16 candidates are structurally clean**, so roughly four
in five histidine combinations on this chassis do not survive refolding whichever positions
they use — a property of the chassis, not of A38H.

Every pKa is PROPKA on the mutant in its bound context, and every combination was measured as a
combination rather than composed from singles — which was necessary: `A35H+A38H` drops A35 from
6.07 to 5.94, out of the window, and `A70H+A72H` drops A70 from 6.44 to 6.23.

| deviation | what was changed, and when |
| --- | --- |
| **Burial floor lowered 15 Å² → 8 Å²** | Changed **after** seeing that `helix_0` had only two switchable positions at the original floor. It then yielded six, four in the ideal window. The floor's purpose — the histidine must sit at the interface for its protonation to modulate binding — is unchanged, but sites admitted between 8 and 15 Å² are more weakly coupled, so "four ideal sites" is a post-hoc number. |
| **`free_4` dropped on a non-reproducing pKa** | Its A16H scanned as 6.58 "ideal" in the first pass and **1.73** on re-measurement. The method is deterministic — `helix_0` and `free_1` reproduce to 0.01 on repeat — so this is not run-to-run noise; the two passes ran against *different structure files* for that design. Treated as disqualifying rather than averaged. |
| **The novelty rule is empirical, and now n = 16 without becoming calibratable** | "High fold similarity fails, moderate passes" is read off Proteinbase verdicts, not from their published scale — which, applied literally at 38% identity, would have failed designs that passed. Sixteen verdicts have not fixed this and will not: **fifteen of the sixteen are the same score (3/4)**, with one 2 and one 4, so there is almost no outcome variance to fit against. Worse, the one controlled pair available — `free1_T73H` and `free1_parent`, one substitution apart with identical TM-score, identical nearest hit and identical identity — scored **4/4 and 3/4**, so the features we measure have *no* resolving power at the only boundary where we do have variance. More verdicts of this kind add nothing. It is used to rank, and it is not calibrated. |
| **The structural `flag` gate was added after seeing the geometry table** | `submission.candidate_problem` originally gated on clash count alone. The first slot-filling geometry table showed two candidates at `epitope_focus 0.00` — bound somewhere else on the target entirely — admitted by a clash count of zero. The flag check was wired in after reading that table, so "one of eight survived" is **not** a pre-registered count. It is not a threshold fitted to the answer: `scripts/submit.py shortlist` has rejected on this same flag since the first campaign, and the omission was in failing to carry it into the newer gate. |
| **`tnf2_13` was added after the first upload** | The twelve-design CSV had already been uploaded and novelty-checked when the eight slot-filling candidates were folded. The design that survived was appended and the whole set re-checked, so `tnf2_13` was selected under gates fixed before it was folded but after the rest of the set was already fixed. |
| **The clash gate was recalibrated after the rounds it judged** | `MAX_CLASHES` was set to **0** for both slot-filling rounds without checking what the designs already submitted score. Seven of the thirteen carry clashes, up to three on `tnf2_04`, so candidates were being held to a standard the incumbents do not meet. Re-gated at that tolerance, `helix0_E38H_E67H_E70H` and `helix0_R35H_E70H_R72H` clear every gate and enter as `tnf2_14` and `tnf2_15`. **This is a threshold changed after seeing which designs it rejected**, so "three of sixteen survived" is not a pre-registered count. It was raised to match the submitted set rather than to admit a particular candidate, and the two it admits were ranked on criteria fixed before they were folded — but the direction of the change was learned from the result. |
| **BindCraft2 runs on loosened design settings** | Rounds 1 and 2 scored zero candidates: trajectories died at the screen stage on binder pLDDT 0.59 against a default floor of 0.60, so the relaxed acceptance filters never saw a candidate. Round 3 applies `target_flexibility 0.5` + `initial_guess` — rungs 1–2 of BindCraft2's own desperation ladder, which cannot self-engage here because it waits 750 trajectories without an acceptance and one trajectory costs ~33 min. **Its documentation is explicit that designs accepted this way are judged against an easier task and are weaker candidates.** Anything this produces is marked as such rather than entering the funnel as an equal. |

## 7. What is not established

- No expression, no SEC, no thermal stability, no SPR/BLI. **No measured binding of any kind.**
  No Tm predictor and no learned solubility model was run — touchstone's
  `ThermostabilityVerifier` needs TemStaPro, which is not installed here, and an uncalibrated
  Tm on de novo miniproteins is not a measurement. What follows is a liability register from
  sequence plus one language-model proxy, not a prediction of behaviour.
- **touchstone's `ExpressionVerifier` returns `trust` for all thirteen, which is not evidence.**
  Its thresholds are library defaults — `max_perplexity` 12.0, `min_solubility` 0.4 — and this
  set runs 6.3–10.5 and 0.51–0.78, so nothing here can fail it. A gate that passes 13/13 at
  untouched defaults has measured nothing; only its *score* is usable, and that score (0.77–0.78
  for `helix_0`, 0.51 for `free_1`) just restates the chassis split below. The solubility half is
  touchstone's own declared placeholder (GRAVY + charge, pending NetSolP/CamSol), so half the
  input is a heuristic and the thresholds were never calibrated against an expression outcome.
- **ESM-2 pseudo-perplexity separates the two chassis and nothing else.** `helix_0` sits at
  6.29–7.26 (mean 6.79, n=11) and `free_1` at 10.10–10.48 (n=2): no overlap, and a gap roughly
  ten times the spread any substitution produces. Lower is more natural, so `free_1` is the
  less protein-like of the two by this proxy — which compounds its 41–47 % mouse retention and
  its near-neutral, hydrophobic surface. Its one advantage, the campaign's only 4/4 novelty,
  is a gate already cleared at 3/4 by everything else. Novelty and naturalness are measured
  against overlapping databases and pull against each other; neither is calibrated to
  expression here.
- **Both liabilities this card previously registered are refuted against measured outcomes.**
  Proteinbase's full export carries `expressed` as an experimental boolean for **4065 designs**,
  which is the first wet-lab calibration this campaign has had, and it contradicts both claims:
  - *The unpaired cysteine.* Eleven of the thirteen carry one cysteine, inherited from position
    34 of the `helix_0` chassis, and this card called the odd count a liability because it
    cannot form a disulphide. Oddness does nothing: designs with an odd count express at
    **86.7 %** and those with an even non-zero count at **86.7 %**, identically. What matters is
    how many, and the effect is small until there are several — within the 60–100 aa band these
    designs occupy, one cysteine costs **0.7 percentage points** (95.5 % against 96.2 % for
    none, n=133 and n=1018). Two or more does cost real yield (88.7 % in the same band, 69.7 %
    above 150 aa), which retrospectively supports dropping `free_4` for carrying three — but
    one cysteine is not the liability this card claimed.
  - *The aggregation axis.* Expression is flat in hydrophobicity across the whole export
    (89.7–93.2 % from GRAVY −2.0 to +2.0, and not monotonic), so `free_1` being near-neutral and
    hydrophobic does not predict that it fails to yield material.

  The feature that does predict expression is **length**, which this card never mentioned:
  96.8 % at 60–80 aa, 95.2 % at 80–100, 91.5 % at 100–150, and 79.6 % above 150. At 74 and 91
  residues both chassis sit in the best bands, against a base rate of **91.6 %**. So the honest
  expectation is that essentially all of these express, and that expression is not where this
  submission is likely to fail. Caveats: the export pools competitions, targets and expression
  systems; these are single-variable rates with one confound controlled, not a model; and none
  of these designs is in that dataset.
- **touchstone's `ExpressionVerifier` returns `trust` for all thirteen, which is not evidence.**
  Its thresholds are library defaults — `max_perplexity` 12.0, `min_solubility` 0.4 — and this
  set runs 6.3–10.5 and 0.51–0.78, so nothing here can fail it. A gate that passes 13/13 at
  untouched defaults has measured nothing; only its *score* is usable, and that score (0.77–0.78
  for `helix_0`, 0.51 for `free_1`) just restates the chassis split below. The solubility half is
  touchstone's own declared placeholder (GRAVY + charge, pending NetSolP/CamSol), so half the
  input is a heuristic and the thresholds were never calibrated against an expression outcome.
- **ESM-2 pseudo-perplexity separates the two chassis and nothing else.** `helix_0` sits at
  6.29–7.26 (mean 6.79, n=11) and `free_1` at 10.10–10.48 (n=2): no overlap, and a gap roughly
  ten times the spread any substitution produces. Lower is more natural, so `free_1` is the
  less protein-like of the two by this proxy — which compounds its 41–47 % mouse retention and
  its near-neutral, hydrophobic surface. Its one advantage, the campaign's only 4/4 novelty,
  is a gate already cleared at 3/4 by everything else. Novelty and naturalness are measured
  against overlapping databases and pull against each other; neither is calibrated to
  expression here.
- **Eleven of the thirteen carry one cysteine, an odd and therefore unpairable count.** It is
  at position 34 of the `helix_0` chassis, so every design on that chassis inherits it rather
  than choosing it. A lone surface cysteine invites disulfide scrambling, covalent dimers and
  heterogeneous material. The campaign already treated this as a liability elsewhere — a
  chassis was dropped partly for carrying three unpaired cysteines — so the eleven ship with a
  known risk that was measured and not acted on. Mutating it to serine was never tested, and
  would need the whole pKa and mouse pipeline re-run to be submittable.
- The two chassis carry **opposite** solubility risks, and neither was measured. `helix_0` is
  strongly acidic and hydrophilic (net −10 to −14 e at pH 7.4, GRAVY −0.37), which is far from
  its pI and argues for good solubility. `free_1` is near-neutral and hydrophobic (−2.2 e,
  GRAVY +0.26) — the combination that aggregates. `tnf2_09` and `tnf2_10` are therefore the
  pair most at risk of never yielding usable material, and they are also the matched control.
- pH selectivity rests on predicted pKa. No folding model represents protonation; PROPKA is a
  static-structure predictor, and the assay is not simulated.
- **The structural ensemble does not explain the pKa spread, and the folding model's own
  confidence cannot be used to screen for it.** Across the eight slot-filling candidates, the
  largest pKa spread per design ranges 0.99–6.06 units while binder pTM varies by 0.002–0.026
  and complex pLDDT by 0.008–0.017 across the same three samples; the correlation between them
  is **−0.46 (n=8)**, i.e. if anything the wrong sign. The extreme case is
  `helix0_R21H_R35H_E67H` position 35, whose pKa reads **3.63, and 9.69**, on structures the
  model scores at binder pTM 0.957 / 0.958 / 0.956. Three structures it considers
  indistinguishable put one histidine anywhere from permanently deprotonated to permanently
  protonated. PROPKA's sensitivity here is to local side-chain placement, not to fold
  uncertainty — so a confidence filter is no substitute for re-scoring every sample, and three
  samples undersample whatever is actually varying. **Every ΔQ floor in this submission is a
  lower bound estimated from three draws, and the true worst case is likely worse.**
- **The structural gate is a coin flip, and every count in the slot-filling funnels is a
  statement about `model_0` rather than about the design.** Clashes and interface area were
  scored on one diffusion sample. Re-scored on all three, **10 of the 16 candidates change
  their clash-free verdict depending on which sample is read**, and interface area moves by a
  factor of two to five on the same design — `helix0_R35H_E67H_E70H` reads 0 clashes over
  658 Å² on `model_0` and 9 clashes over 1203 Å² on `model_2`. The parent chassis is clash-free
  on its own sample, and a single histidine is enough to introduce a clash on three of the four
  one-histidine designs, so the gate is penalising the substitution the design exists to make.
  Numbers in `docs/clashes-per-sample.txt`.

  This is the same mistake as the pKa, made a third time. A quantity was measured on one
  arbitrary sample and treated as a property of the design, after that exact error had already
  been found and fixed for pKa (`switch.SampledHis`) and for interface confidence (mean rather
  than best-of-n). **It is not fixed here**: the structural verdicts in this card stand as
  `model_0` measurements, and the honest reading of "six of eight failed on structure" is "six
  of eight failed on the one sample that was scored".
- **The three samples are not energetically equivalent, and the ΔQ floor is usually set by the
  worse one.** Rosetta REF15 on the raw coordinates separates them by 347–2855 REU across the
  eight slot-filling candidates — the folding model rates the same triples at pTM 0.956–0.958.
  Taking each histidine at its least favourable sample, as `delta_q_floor` does, draws on the
  best-scoring structure for only **9 of 24 histidines**; the other 15 come from structures
  REF15 ranks worse, in one case by 2855 REU. The headline example: `helix0_R21H_R35H_E67H`
  His35 reads pKa **9.69 on the structure scored 731 REU and 3.63 on the one scored 64.8**.
  So the ΔQ floors here are systematically *too* pessimistic, by 0.02–0.59 per design.
  **It did not change which design was selected.** Re-ranking the eight on their best-scoring
  structure promotes three others above `helix0_R35H_E67H_E70H` — and all three had already
  been dropped by the structural gate, for clashes or an interface off the epitope. Of the two
  that cleared that gate, `R35H_E67H_E70H` leads under both rules (1.35 vs 0.92 worst-sample,
  1.60 vs 0.94 best-structure). Checked rather than assumed; raw REF15 on predicted coordinates
  is dominated by unrelaxed clashes, so the absolute numbers are ranking-only, but a 2855 REU
  gap is far outside what a relax would absorb.
- Mouse cross-reactivity rests on a *predicted* mouse structure and conserved-epitope coverage.
- Three mechanistic hypotheses for the pH switch (buried-His desolvation, cation repulsion,
  carboxylate proximity) were each **refuted** by the pKa data; no validated mechanism replaces
  them. Positions come from a scan, not a rule.

## 8. Stop rules

Stop and re-think if: no construct reaches fold-novelty ≥3; the His scan yields no pKa in
6.3–7.2 at an interface position; or cross-model interface confidence stays below ~0.3 for
every design, which would mean the ranking cannot enrich anything.
