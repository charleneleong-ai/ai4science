# Conditional TNF-α binder — supporting information

Anthropic × Adaptyv Challenge 02, Track 3. Twenty single-chain de novo designs, lengths 74–91 aa.

`molecule_class` is `single_chain`. Worth noting that the documentation and the uploader differ
on this: FAQ 5 and the challenge page both enumerate the permitted values as `protein`,
`nanobody`, `scfv`, `fab_kappa`, `fab_lambda`, while the downloadable submission template labels
its single-chain example `single_chain`. `single_chain` is what the uploader accepted for the
first submission, so it is what is used here.

Everything here is computational. **No binding of any kind has been measured** — no expression,
no SEC, no thermal stability, no SPR or BLI. Each number below carries the sampling depth it was
measured at, because on this campaign the depth changed the answer.

## Files

| file | what it is |
| --- | --- |
| `metrics.csv` | 20 designs × 24 columns, one row per submitted design, in submission order |
| `methodology/method.md` | epitope selection, the pH-switch mechanism, and why the obvious mechanism fails |
| `methodology/experiment-card.md` | frozen target, exact tool versions, the funnel with denominators, declared post-hoc filter changes |
| `methodology/sampling-depth.md` | the main methodological result: what three diffusion samples could not tell us |
| `methodology/architecture.mmd` | the pipeline as built, with per-stage survivor counts |
| `methodology/report.pdf` | the written report |
| `structures/rank{NN}_{name}.{cif,pdb}` | predicted complex per design, named by submission rank |
| `structures/manifest.csv` | rank → design → file mapping |

Code and full history, public: <https://github.com/charleneleong-ai/ai4science/tree/main/tnf-binder>

## Ranking method

Submission order follows the challenge's own objective order — pH-selectivity first, mouse
cross-reactivity second, affinity third — with each term taken at the deepest sampling available:

1. **Number of switching histidines.** A histidine counts as switching if its predicted pKa falls
   in the 6.0–7.6 window on ≥90% of **20** diffusion samples. Within an equal count, a design
   measured at 20 samples ranks above one still provisional on 3, so that an unmeasured design is
   neither credited nor penalised as though it had been tested.
2. **Mouse retention**, the ratio of symmetrised binder-to-target ipTM means, mouse over human.
3. **Human binder-to-target ipTM**, mean over samples rather than best-of-n.

The 0.90 rate threshold is **post hoc** — chosen after seeing the distribution. The full
sensitivity curve is reported rather than the single cut: at a rate of 1.00 it selects 0 designs
with ≥2 switching histidines, 0.95 selects 2, **0.90 selects 5**, 0.85 selects 7, 0.80 selects 9,
0.75 selects 10.

## The result that matters most

The pre-registered criterion was "in the pH window on every sample", measured on **three**
samples. Re-measuring 25 designs at **20** samples (500 structures, one PROPKA run per histidine
per structure) changed three things:

- **Two of the four positions the campaign was built on are not switches.** Pooled in-window
  rates: His67 0.94, His72 0.84, His70 0.83, but His21 **0.56**, His38 0.56 and His35 **0.43**,
  with median pKa spreads of 4.28 and 5.70 pH units between structures the folding model scores
  equally. A predicted pKa with that spread is a property of the sample, not of the design.
- **"In window on every sample" is a maximum**, so it tightens as samples are added and cannot be
  compared across sample counts. At n=3 most designs passed a ≥2-histidine gate; at n=20, **0 of
  25** do. The comparable quantity is the in-window *rate*, which is what is used here.
- **The rate is design-specific, not positional.** His70 and His72 hold up individually but
  degrade to 0.65 and 0.60 inside `E67H_E70H_R72H`, so switch counts cannot be composed from
  single-substitution scans.

Mouse retention also moved by up to 0.25 **in both directions** between 3 and 20 samples, which
reversed which designs ranked best. Because the submission order turned on it, the three
top-tier designs that still had only 3-sample mouse data were re-folded at 20 samples on
2026-10-17. All three 3-sample values were optimistic — 0.874 → **0.802**, 0.824 → **0.735**,
0.853 → **0.860** — which reordered ranks 3 to 5. **All five top-ranked rows now rest on 20
samples.** Six rows further down still carry 3-sample retention and are flagged per row in
`metrics.csv`; two of them (`tnf2_06`, `tnf2_07`) were never folded at 20 samples at all.

## Funnel

Failed generations are in the denominator. Stage counts are in
`methodology/experiment-card.md` §4; the final stage is summarised here.

| stage | n |
| --- | --- |
| BoltzGen designs generated | 100 |
| campaign constructs (parents + histidine variants) | 57 |
| chassis passing both the fold-novelty and composition gates | **2** of 16 folded parents |
| designs folded at 20 samples | 25 human / 28 mouse (500 + 560 structures) |
| designs with ≥2 switching histidines at rate ≥0.90 | **5** |
| designs with ≥1 | 18 |
| submitted | 20 |

The concentration on two chassis is a result of the funnel, not a preference: of 16 folded
parents, only `screen_pocket_helix_0` and `screen_pocket_free_1` passed both gates. 18 of the 20
rows are a measured histidine series on the first chassis; pairwise identity among them is
0.934–0.989, so they are a substitution series rather than independent designs.

## Stack

| role | tool | version / checkpoint |
| --- | --- | --- |
| design (sampling) | BoltzGen | 0.1.4, `boltzgen1_adherence` + `boltzgen1_diverse` |
| design (optimisation) | BindCraft | 1.5.3, AF2 `params_2022-12-06` |
| design (dual-target) | BindCraft2 | 1.0.3, jax 0.11.2 |
| co-folding / validation | Boltz-2 | 2.2.1, `boltz2_conf.ckpt` |
| side-chain building | PyRosetta | 2026.29 |
| pKa prediction | PROPKA | 3.5.1 |
| sequence novelty | MMseqs2 | vs Swiss-Prot + PDB seqres, snapshot 2026-10-10 |
| fold novelty | FoldSeek + TM-align | CATH50, snapshot 2026-10-11 |

Boltz-2's affinity head is not used anywhere: it requires the binder chain to be a small-molecule
ligand and caps at 128 atoms, against ~450 heavy atoms for these designs.

## Per-design caveats

- `tnf2_18` — its campaign `model_0` read 232 Å² with a tiny-interface flag; the 20-sample median
  is 825 Å² and never below 535, so the single model was the artefact. The structure and geometry
  supplied here come from the deep fold (941 Å², 4 clashes), not from that model.
- `tnf2_20` — the uploaded structure is the `sub2` fold while its single-model geometry came from
  the `a38` scan; the 20-sample geometry in `metrics.csv` supersedes both.
- `tnf2_06`, `tnf2_07` — never folded at 20 samples. Their switch counts are provisional on 3
  samples and they are ranked below every design measured at 20 with the same count. Blank
  20-sample geometry columns mean not measured, not zero.
- `tnf2_09`, `tnf2_10` — the free_1 chassis. `tnf2_09`'s His73 has a mean pKa of 5.94, below the
  window floor, in-window on 11 of 20 samples, so it is **not** a switch; `tnf2_10` carries no
  histidine at all. Both are ranked last on the primary objective. At 20 samples their human
  binder-to-target ipTM is 0.434 and 0.424 — roughly 2.4× the helix_0 series, not the ~4× the
  3-sample numbers suggested — and their mouse retention of 0.466 and 0.482 is the worst in the
  set. They are included as the only second chassis the funnel produced.

## What is not established

- No measured binding, affinity, expression or stability. Everything is prediction.
- pH selectivity rests on PROPKA over predicted structures. No folding model represents
  protonation, and the pH 7.4 / 6.0 assay is not simulated.
- Mouse cross-reactivity rests on a predicted mouse trimer and conserved-epitope coverage, not on
  an experimental mouse structure.
- Three mechanistic hypotheses for the switch — buried-histidine desolvation, cation repulsion,
  and carboxylate proximity — were each refuted by the pKa data. Positions come from a brute-force
  scan plus a physics solver, not from a rule, and no validated mechanism replaces them.
- Interface confidence is not a point estimate. On identical structures AF2 scored a BindCraft
  design 0.92 where Boltz-2 scored 0.21, so no single model's number is reported as truth, and no
  model scores designs it also optimised.
- The filter set was not fixed before results were seen. The composition gate, the fold gate and a
  lowered burial floor were each added in response to an observed failure; all are declared in
  `methodology/experiment-card.md` §6 and §6b. The funnel above is therefore not a clean
  pre-registered statistic.
