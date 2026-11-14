# Spectrum-to-structure scorer options (2026-11-14)

Goal: a better or complementary scorer to z-blend with sanpier v1 + popwenman v80 within ~4 weeks.
Constraints: Kaggle T4, internet off (weights shipped as a dataset), A100 for training. Host rules and forum LB
evidence are in [kaggle_discussion.md](kaggle_discussion.md). Numbers are author-reported unless labelled "audit"
([MassSpecGym in the Wild, arXiv 2606.19624](https://arxiv.org/abs/2606.19624)). MSG = MassSpecGym. T4 costs are our
estimates unless a timing is quoted. Nothing was downloaded or run.

## What is new since 2026-10-30

1. **MassSpecGym is cleared as training data** ("MassSpecGym and MIST models are fine", host, 11-11,
   [757307](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/757307)). This moves
   every MSG-only checkpoint (FLARE, JESTR, MIST `msg_*`, GLACIER/ICEBERG-MSG) from "unclear" to allowed, provided its
   training data is MSG-only and the licence permits redistribution.
2. **DreaMS official checkpoints are allowed** (host, 11-05,
   [754610](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754610)); final
   submissions must be "reproducible end-to-end ... Yes".
3. **FIORA-OS v1.0 is the repo default and is trained on MSnLib only.** The [BAMeScience/fiora](https://github.com/BAMeScience/fiora)
   README (MIT, last push 2026-04-28) says an open-source model is selected by default and `fiora-train` reproduces it
   from MSnLib v7 on Zenodo. The Nat Commun paper ([PMC11889238](https://pmc.ncbi.nlm.nih.gov/articles/PMC11889238/))
   says the *default published* model is NIST'17 + MS-DIAL (ineligible) and FIORA-OS is 80% of MSnLib v1.0 (177k
   spectra, 16k compounds).
4. **MS2DeepScore dual-mode v5** ([Zenodo 13897744/latest](https://zenodo.org/records/13897744/latest), 2025-12-05,
   CC BY 4.0, 416 MB) adds MSnLib spectra and extra cleaning; training data at Zenodo 16882111.
5. **No new primary-source scorer paper or HF checkpoint found** since 10-30. HF model search for MassSpecGym shows
   nothing newer than 2026-07-30; MassSpecGym repo last push 2026-05-08; JESTR1 2026-05-08; DreaMS 2026-06-05. The arXiv
   API returned HTTP 429 twice, so the arXiv delta is incomplete (see "Not verified").

## Candidates

### 1. FLARE (Hassoun lab) — fine-grained peak↔atom contrastive scorer
- Primary: [PMC12873900](https://pmc.ncbi.nlm.nih.gov/articles/PMC12873900/) (bioRxiv 2026.01.27.702086). Weights
  `pretrained_models/flare.ckpt` in [HF Space HassounLab/FLARE](https://huggingface.co/spaces/HassounLab/FLARE), MIT,
  trained on MSG only (per 10-17 sweep).
- Numbers: MSG mass track hit@1/5/20 43.15/75.59/92.89 (author); formula track 22.66/50.00/75.15 (audit confirms).
  Best audited *contrastive* scorer; well above JESTR/MVP.
- Rules: MSG now allowed, MIT → allowed (inferred from 757307; no FLARE-specific ruling).
- Needs subformula-annotated peaks (MIST-style), so a formula per candidate: we have it from the candidate SMILES.
- T4: small GNN + transformer; 400 × ~200 candidates = 80k graph encodes, est. minutes.
- timsTOF fit: none (MSG is Orbitrap/QTOF); retrain/fine-tune on train.parquet on the A100 is possible (code public).
- Family: different from the kit scorer (fine-grained alignment vs global embedding) → best complementarity bet.

### 2. MIST `msg_aug` fingerprint scorer (coleygroup)
- [coleygroup/MIST](https://github.com/coleygroup/MIST) (MIT, training code), [Zenodo 23303289](https://zenodo.org/records/23303289).
  Host-allowed by name (751373, 757307).
- Numbers: no retrieval number for `msg_aug` in the record; audit formula hit@1 for MIST vMSG 9.57 (older checkpoint).
  Our probe: mean FP Tanimoto 0.385 on 50 Enveda-180 timsTOF spectra vs 0.261 for the old FRIGID MIST (failed our
  0.40 de-novo gate, but that gate was for generation, not ranking).
- Scorer use: cosine / log-likelihood between predicted FP and candidate FP, z-blended. Not tested as a ranker.
- T4: ~290 MB; ~0.6 s/spectrum (batch 1 because of the attention-mask bug, audit §6.1.1); 400 molecules trivial.
- timsTOF fit: none; can fine-tune on train.parquet with MIT code.
- Note: the public lineage already has an FPNet; independent FP engines added +0.001 in 749251 and 7 FPNet swaps lost
  (kaggle_discussion.md). Blend only, never replace.

### 3. DreaMS fine-tuned as a retrieval head
- [pluskal-lab/DreaMS](https://github.com/pluskal-lab/DreaMS) (MIT); weights [HF roman-bushuiev/DreaMS](https://huggingface.co/roman-bushuiev/DreaMS),
  Zenodo 10997887 (CC BY 4.0). Host: official checkpoints allowed (754610); fine-tune of `ssl_model.ckpt` asked but
  unanswered.
- Numbers: no primary MSG retrieval number for a fine-tuned DreaMS. Frozen DreaMS + frozen molecule encoder alignment
  caps at MRR 0.227 on decontaminated MSG v1.5 ([arXiv 2610.37979](https://arxiv.org/abs/2610.37979)); audit corrected
  DreaMS+ChemBERTa to 6.70 hit@1. So the value must come from fine-tuning the encoder on our 2.4M spectra.
- Pretraining corpus GeMS is mined from MassIVE/GNPS and includes many instrument types, plausibly some timsTOF
  (not verified).
- T4: ~116M-param transformer (not verified), 400 × few spectra = trivial; candidates via a molecule tower.
- Effort: highest of the six (train a two-tower retrieval head on A100, ~3-5 days); also the most different
  representation from the kit scorer.

### 4. FIORA-OS forward model → cosine to the query
- See item 3 above. MIT code; OS weights in-repo; trained on MSnLib (open). Host-allowed by name (751373).
- Numbers: no top-k identification numbers in the main text; median sqrt-cosine 0.81 on its own test split vs ICEBERG
  0.72 and CFM-ID 0.67; on CASMI 2022 0.29 vs CFM-ID 0.38 and ICEBERG 0.36 (author; default NIST model, not OS).
- Speed: 0.03 s/spectrum on A100, 0.92 s on CPU (paper Table 2). T4 est. ~0.1 s → 80k candidate spectra ~2 h;
  restrict to same-formula top-25 (~10k) → ~15 min.
- CE input is continuous eV; test CE is nominal (753419). Instrument type is a covariate (no timsTOF class known).
- Forum: a FIORA-OS negative-mode blend at 0.25 lost −0.014 offline panel MRR (zz3t3rn41 lineage). Overlaps GLACIER
  (already in the pipeline). Low priority.

### 5. MS2DeepScore 2.0 (dual-mode v5) — spectrum↔spectrum similarity
- [matchms/ms2deepscore](https://github.com/matchms/ms2deepscore) (Apache-2.0, push 2026-09-15); weights
  [Zenodo 13897744 v5](https://zenodo.org/records/13897744/latest) CC BY 4.0. Host-allowed by name (751373).
- Not a structure scorer: it predicts structural similarity between two *spectra*. Uses here:
  (a) analogue-library scorer vs the extlib / MassBank references (class 1/2), (b) query vs FIORA/GLACIER-predicted
  candidate spectra, (c) test↔test edges for transductive propagation.
- Numbers: no MSG retrieval number (it is not evaluated that way). Cross-mode preprint
  [bioRxiv 2024.03.25.586580](https://www.biorxiv.org/content/10.1101/2024.03.25.586580) — numbers not pulled.
- T4: embedding a few thousand spectra is seconds; 416 MB weights.
- timsTOF fit: trained on GNPS/MassBank/MSnLib (mostly QTOF/Orbitrap); the 756504 finding that timsTOF↔other
  instruments sit at ~0.90 entropy similarity suggests a learned similarity may transfer better than cosine.

### 6. JESTR (Hassoun lab) — global contrastive
- [HassounLab/JESTR1](https://github.com/HassounLab/JESTR1) (MIT); MSG weights in `data/MassSpecGym/`
  (`pretrained_spec_enc_model_*`, `pretrained_mol_enc_model_*`, `pretrained_inter_model_*`). MSG → now allowed.
- Numbers: audit formula 11.82/33.48/61.46; mass 15.13/36.75/60.32 (FLARE table).
- Forum: RRF with JESTR on the 0.440 base gave 0.418 (−0.022; 2026-10-17 sweep). Cheapest drop-in but already lost once
  (as RRF, not as a z-blend at ranks 2-25).

## Excluded or low priority

| Model | Reason |
|---|---|
| MVP ([HassounLab/MVP](https://github.com/HassounLab/MVP)) | no LICENSE file (redistribution rule); weaker than FLARE (audit formula 13.96 hit@1) |
| CSU-MS2 ([tingxiecsu/CSU-MS2](https://github.com/tingxiecsu/CSU-MS2)) | no LICENSE; weight location/training data not stated in README (unverified); no MSG numbers |
| MassFormer | published models trained on NIST (ineligible); would need retraining (not checked further) |
| FraGNNet | MSG mass/formula hit@1 46.64/31.93 (GLACIER table), but no official weights found (10-17 sweep) |
| ICEBERG / SCARF | ICEBERG already in the pipeline; SCARF superseded, no weights |
| MS-BART, DiffMS, MADGEN | de novo; audit top-1 ≤ 2.3% |
| GLMR, SpecBridge | shortcut-inflated per audit; no code/weights |
| Separability recipe (2610.37979) | best decontaminated MRR 0.566, but no weights; 4×A100 training plus simulating pools; possibly the strongest option if time allowed |
| SIRIUS / CSI:FingerID | banned (751373) |
| Kit retrains on Kaggle (bgunshui v2-v8, harshitkmr b160w, studyexchange real-adduct, youyouwou v9np, czhuangsjtu phys) | same family as what we blend; forum shows ±0.003 at 50/50 and losses when weighted up |

## Not verified

- FLARE's licence and data provenance were taken from the 10-17 sweep, not re-read today; no host ruling names FLARE.
- Whether FLARE's mass-track number used GT or predicted formulas.
- DreaMS parameter count, GeMS timsTOF share, and the host answer on fine-tuning `ssl_model.ckpt`.
- MIST `msg_aug` as a candidate ranker (only FP Tanimoto measured).
- MS2DeepScore v5 benchmark numbers; FIORA-OS (vs default) identification numbers; all T4 runtimes.
- The arXiv delta since 10-30 (API rate-limited); the current leaderboard top (MCP unauthenticated).
