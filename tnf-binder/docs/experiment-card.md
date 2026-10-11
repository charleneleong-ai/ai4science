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

## 7. What is not established

- No expression, no SEC, no thermal stability, no SPR/BLI. **No measured binding of any kind.**
- pH selectivity rests on predicted pKa. No folding model represents protonation; PROPKA is a
  static-structure predictor, and the assay is not simulated.
- Mouse cross-reactivity rests on a *predicted* mouse structure and conserved-epitope coverage.
- Three mechanistic hypotheses for the pH switch (buried-His desolvation, cation repulsion,
  carboxylate proximity) were each **refuted** by the pKa data; no validated mechanism replaces
  them. Positions come from a scan, not a rule.

## 8. Stop rules

Stop and re-think if: no construct reaches fold-novelty ≥3; the His scan yields no pKa in
6.3–7.2 at an interface position; or cross-model interface confidence stays below ~0.3 for
every design, which would mean the ranking cannot enrich anything.
