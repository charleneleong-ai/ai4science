# MS/MS -> structure methods sweep (2026-10-17)

Raw notes for CASMI 2026 (Enveda, Kaggle). All numbers are **vendor-reported** unless marked
"audit" (= independent re-evaluation in MassSpecGym in the Wild). MSG = MassSpecGym.
"(unverified)" = could not confirm from a primary source. T4 cost lines are our own estimates
unless a timing is quoted.

## 0. Read this first: the MSG audit

- **MassSpecGym in the Wild** (Liu, Bushuiev, ..., Hassoun, Pluskal, Coley; arXiv 2606.19624, Jun 2026)
  https://arxiv.org/abs/2606.19624 , HTML https://arxiv.org/html/2606.19624
  - Evaluation issues in >=17 of 26 papers reporting MSG results in year one: data leakage, shortcut
    learning, implementation bugs / metric divergence. Releases MSG v1.5 (https://github.com/pluskal-lab/MassSpecGym, MIT).
  - SMILES-canonicalization shortcut (non-canonical ground truth vs canonical decoys):
    DreaMS+ChemBERTa alignment 82.41 -> 6.70 hit@1; ChemFormer two-stage 63.89 -> 7.60;
    ChemBERTa format classifier 99.66 -> 8.2. (These designs correspond to SpecBridge and GLMR —
    our inference from architecture; the audit excerpt read did not name papers.)
  - MIST batched-inference bug: boolean padding mask added to attention (+1 to pads) inflates FP
    predictions (§6.1.1). MIST vDiffMS checkpoint not reproducible. MIST+MolForge 31.75% -> 10.73% top-1 corrected.
  - Formula leakage, mass track (Table 2, hit@1/5/20): random 0.43/2.01/7.99; random **with GT formula
    12.98/39.09/71.33**; MIST + MIST-CF top-1 formula 26.64/38.47/54.90; MIST + GT formula
    **51.45/66.25/80.58**. => correct formula roughly doubles hit@1.
  - Spectrum-blind PubChem default ordering on formula track: 49.98/54.13/61.04 hit@1/5/20 —
    popularity is a huge prior (supports our popularity prior; also means decoys must be popularity-matched).
  - Formula-track retrieval (Table 3, hit@1/5/20): Random 3.06/11.35/27.74; MIST vMSG 9.57/22.11/41.12;
    JESTR 11.82/33.48/61.46; MVP 13.96/36.88/68.12; FLARE 22.66/50.00/75.15; ICEBERG 36.18/60.70/78.83;
    DiffMS (generative retrieval) 32.29/54.36/75.47; MIST vDiffMS 46.42/54.74/66.27 (buggy ckpt);
    FRIGID 45.27/58.30/69.46; **MIST vFRIGID 53.76/65.32/74.71** (not stated whether re-run). No MRR reported.
  - Valid de novo leaderboard (Table 9, top-1 / top-10 acc, top-1/top-10 Tanimoto): MSNovelist 0/0;
    Spec2Mol 0/0; MADGEN 1.31/1.54; MS-BART 1.07/1.11; DiffMS 2.30/4.25 (0.28/0.39);
    FOAM 1.50/10.28 (0.35/0.53); MIST+MolForge corrected 10.73/14.48 (0.37/0.41); **FRIGID 18.29/22.00 (0.43/0.47)**.
  - Recommendations: uniform RDKit canonicalization; spectrum-masking control; property-matched decoys;
    keep formula-known vs mass tracks separate; simulators trained on train split only;
    14-char InChIKey matching, 2048-bit ECFP4, bootstrap CIs.

## 1. Strongest MSG retrieval as of 2026-10

- Highest number found: **GLACIER** (forward model) mass track hit@1/5/20 **69.95/86.47/93.52**, formula
  track 49.88/78.13/90.75 (own paper, Table 1; up to 256 candidates/spectrum, official MCES split).
  https://arxiv.org/html/2606.29161 — vendor-reported, not in the audit (posted after).
- Caveat: the same Table 1 reports ICEBERG 2.0 at 44.35 formula hit@1, while the audit measured ICEBERG
  at 36.18. Forward-model retrieval swings with collision-energy imputation and candidate caps.
- Highest audit-table row: MIST vFRIGID 53.76 formula-track hit@1.
- Inflated / discard: SpecBridge 84.7 R@1, GLMR 64-68 hit@1 (shortcut-prone designs).

## 2. Retrieval / embedding methods

### MIST (Goldman et al., Nat Mach Intell 2023)
- https://github.com/samgoldman97/mist (MIT). Spectrum -> fingerprint via subformula-peak transformer.
- Numbers: MSG mass track 14.64/34.87/59.15 (FLARE table, PMC12873900); formula track 9.57/22.11/41.12 (audit);
  MIST+MIST-CF top-5 formula element-wise max 28.23/40.77/56.22 (audit Table 8).
- Weights: quickstart download script; training data of bundled weights not stated; csi2022 data includes NIST
  (not redistributed) -> eligibility unverified; retrain on our 2.5M spectra. MIST vFRIGID encoders on
  Zenodo 19685145 are CC BY-NC 4.0 (see FRIGID).
- Bug: batched attention mask (audit §6.1.1) — run batch=1 or patch.
- T4: 98.6M params, ~0.64 s/instance (GLMR appendix, https://arxiv.org/html/2511.06259v1). 400 queries trivial;
  candidate scoring = FP cosine, free.
- Plug-in: stronger FPNet replacement / second FP feature (cosine, per-bit log-likelihood) for the ranker.

### MIST-CF (Goldman et al. 2023)
- https://github.com/samgoldman97/mist-cf (MIT). Formula + adduct annotation; public model trained on NPLIB1/CANOPUS,
  NIST20 models on request. Top-1 formula accuracy not pulled from primary source (unverified).
- Plug-in: formula posterior feature / pool narrowing. Biggest lever per audit Table 2.

### JESTR (Kalia et al. 2025)
- https://github.com/HassounLab/JESTR1 (MIT). Weights for NPLIB1 and MSG in repo `data/`; no NIST.
- Numbers: mass 15.13/36.75/60.32; formula 11.85/32.95/61.46 (FLARE table); GLMR paper cites MRR 40.4 (mass) / 33.3 (formula).
- T4: 19.9M params, ~0.57 s/instance. Plug-in: contrastive score feature; superseded by FLARE.

### MVP (multiview; Chen & Hassoun 2025)
- https://github.com/HassounLab/MVP. README mentions MSG-trained model; no weights URL / LICENSE found (unverified).
- Numbers: mass 26.37/58.90/86.88; formula 11.10/31.14/61.99 (FLARE table); audit formula 13.96/36.88/68.12.
- Plug-in: ranker feature; superseded by FLARE.

### FLARE (Chen & Hassoun; bioRxiv 2026.01.27.702086)
- https://www.biorxiv.org/content/10.64898/2026.01.27.702086.full.pdf , https://pmc.ncbi.nlm.nih.gov/articles/PMC12873900/
- FILIP-style fine-grained contrastive alignment of subformula peak tokens vs GNN atom nodes.
- Numbers: **mass 43.15/75.59/92.89**, MCES@1 7.84; formula 22.66/50.00/75.15 (audit confirms formula row). No MRR.
- Caveat: inputs are MIST-style subformula peaks (need a formula). Whether mass-track number uses GT formula is
  unverified; repo has an msbuddy-predicted-formula recipe (`run_scripts/msgym_cand_by_mass_msbuddy_pred_form`) — check it.
- Weights: MIT, `pretrained_models/flare.ckpt` in HF Space https://huggingface.co/spaces/HassounLab/FLARE ; trained on MSG only (open data).
- T4: small GNN+transformer; 400 x 5k = 2M candidate graph encodes, est. tens of minutes (unverified).
- Plug-in: ranker feature (fine-grained score) or extra fusion engine.

### SpecBridge (arXiv 2601.17204, Hassoun lab)
- https://arxiv.org/html/2601.17204 . DreaMS encoder (last 2 blocks tuned) + residual mapper onto frozen ChemBERTa-2.
- Numbers: MSG R@1 84.7 (track unclear) — almost certainly the canonicalization shortcut (audit: 82.41 -> 6.70).
  Spectraverse w/ PubChem formula pools (~1494 cands) R@1 36.6; MSnLib 53.4 (possibly cleaner, unverified).
- Code https://github.com/HassounLab/SpecBridge returned 404; licence unknown.
- T4: very cheap. Plug-in: feature only with strict canonicalization; low priority.

### GLMR (arXiv 2511.06259, AAAI 2026)
- https://arxiv.org/html/2511.06259v1 . ChemFormer contrastive pre-retrieval + generative rerank.
- Numbers: mass 64.17/72.96/78.78, MRR 67.82; formula 68.48/78.09/84.22, MRR 72.47 — likely shortcut-inflated
  (audit ChemFormer two-stage 63.89 -> 7.60). No code/weights found. Skip.

### DreaMS (Bushuiev et al., Nat Biotech 2025)
- https://github.com/pluskal-lab/DreaMS (MIT); weights https://huggingface.co/roman-bushuiev/DreaMS (MIT), Zenodo 10997887.
- Self-supervised on GeMS (GNPS/MassIVE mined; no NIST). 1024-d embedding. No primary-source MSG retrieval number found;
  frozen-embedding retrieval is weak (audit corrected DreaMS+ChemBERTa 6.70 hit@1).
- Plug-in: we already have embeddings; best use = init for fine-tuned spectrum->FP head on our 2.5M spectra.

### CSU-MS2 (Xie et al., Anal Chem 2025, doi 10.1021/acs.analchem.5c01594)
- Contrastive on large simulated corpus. No official repo found; Zenodo 10.5281/zenodo.20794027 seen only on a
  third-party page (unverified). No MSG numbers.

### MSBERT (Anal Chem 2024, doi 10.1021/acs.analchem.4c02426)
- Masked + contrastive spectral embedding for library matching; data Zenodo 13722644. No repo/licence/MSG numbers found. Low priority.

### LSM-MS2 (Matterworks, arXiv 2510.26715)
- https://arxiv.org/html/2510.26715v1 . Claims +30% isomer accuracy, +42% IDs in complex samples; no MSG numbers;
  apparently commercial, no open weights (unverified). Skip.

### "Small molecule retrieval from MS/MS: what are we optimizing for?" (arXiv 2602.16507)
- https://arxiv.org/abs/2602.16507 . Theory: FP accuracy vs retrieval trade-off; loss / FP choice should depend on
  candidate-set similarity structure. Suggests retrieval-aware (listwise over candidate pool) FPNet loss.

### Selective prediction for retrieval (arXiv 2603.10950)
- https://arxiv.org/pdf/2603.10950 — confidence estimation; not read for numbers.

## 3. Forward simulators / rescoring

### GLACIER (Coley lab, arXiv 2606.29161, Jun 2026)
- https://arxiv.org/abs/2606.29161 . Single-stage transformer, fragment-subgraph detection.
- MSG Table 1 (mass / formula hit@1/5/20): GLACIER+contrastive FT 69.95/86.47/93.52 / 49.88/78.13/90.75;
  GLACIER no FT 69.69/86.58/93.34 / 48.49/76.56/90.60; ICEBERG 2.0 63.95/83.83/92.56 / 44.35/73.57/89.45;
  MARASON - / 34.03/64.04/85.39; MolSpecFlow 55.32/79.98/89.54 / -; FraGNNet 46.64/72.56/83.58 / 31.93/63.20/82.70;
  NEIMS 8.44/21.43/38.57 / 7.62/22.70/44.12.
- NIST'20 top-1 52.5% (FT) / 38.5% (no FT).
- Speed: 368.7 spectra/s/GPU (2x A5000) vs ICEBERG 2.0 46.4.
- Weights: ms-pred README (MIT) links MSG-trained GLACIER checkpoint (Dropbox); NIST weights only to licence holders.
  https://github.com/coleygroup/ms-pred
- T4: est. 100-150 spectra/s -> 2M preds (400 x 5000 x 1 CE) ~4-5 h; 1000 cands < 1 h.
- Plug-in: already in pipeline. Check (i) contrastive-FT checkpoint, (ii) multi-CE scoring (max over CEs).

### ICEBERG 2.0 / 2.1 (Wang, Manjrekar, ..., Coley; bioRxiv 2025.05.28.656653)
- https://www.biorxiv.org/content/10.1101/2025.05.28.656653v1 . Adds CE and polarity.
- NIST'20 [M+H]+ top-1 40%, ~92% top-10. MSG: see GLACIER table; audit 36.18 formula hit@1.
- Weights: MSG-trained ICEBERG 2.1 (msg_all, msg_simulation) via ms-pred README (MIT code).
- Eligibility flag: FOAM notes MSG ICEBERG had missing CEs imputed by a NIST'20-trained ICEBERG
  (https://arxiv.org/html/2602.07709) — check msg_all lineage (unverified).
- T4: est. ~15 spectra/s -> rescore shortlist only (top-100-300). Plug-in: existing rescoring feature; second opinion to GLACIER.

### MARASON (ICML 2025, arXiv 2502.17874)
- https://arxiv.org/abs/2502.17874 . ICEBERG + retrieval of reference spectra + neural graph matching.
- NIST top-1 28% vs 19% w/o RAG; MSG formula hit@1 34.03 (GLACIER table). Code in ms-pred, no checkpoint. Low priority.

### SCARF (Goldman et al., arXiv 2303.06470)
- https://arxiv.org/abs/2303.06470 . Formula prefix-tree predictor (ms-pred). Numbers not retrieved (unverified).
  No weights linked. Superseded.

### FraGNNet (Young, ..., Röst; arXiv 2404.02360, TMLR 2025)
- https://arxiv.org/abs/2404.02360 . MSG 46.64 mass / 31.93 formula hit@1 (GLACIER table). No official weights found (unverified). Skip.

### ESP (Li & Hassoun, Bioinformatics 2024, btae490)
- https://pubmed.ncbi.nlm.nih.gov/39180771/ , https://arxiv.org/abs/2203.13783 . MLP+GNN ensemble; NIST20-trained
  (ineligible). Avg-rank gains 23.7% / 37.2% over MLP / GNN. Skip.

### FOAM (Coley group, arXiv 2602.07709)
- https://arxiv.org/html/2602.07709 . Formula-constrained graph GA, ICEBERG oracle, seeds = formula-matched PubChem + DiffMS samples.
- NIST'20 top-1 11%, top-10 31%, true structure encountered 68.1%. MSG (1500 oracle calls) top-1 1.50%, top-10 10.28%, top-10 Tanimoto 0.53.
- ~3 min/spectrum on H100/H200. Code https://github.com/coleygroup/foam (MIT per fork B; paper licence CC BY-NC-SA per fork C);
  MSG weights/seed caches Zenodo 18502041 (CC BY 4.0). NIST lineage risk via CE imputation.
- T4: full 400 x 1500 calls ~60-100 h (infeasible); FOAM-lite (~200 calls seeded with our top-50, GLACIER oracle) est. 2-4 h.
- Plug-in: candidate generator (analog mutation around top hits).

### FIDDLE (Hong et al., Nat Commun 2025; bioRxiv 2024.11.25.625316)
- https://www.biorxiv.org/content/10.1101/2024.11.25.625316 , https://msfiddle.readthedocs.io/ . TCN formula predictor + FDR refinement;
  Orbitrap / Q-TOF; optional BUDDY/SIRIUS. Pretrained download exists; training data / licence not found (NIST plausible, unverified).
- Related formula-tool comparison: arXiv 2601.00941. Plug-in: formula prior/feature.

### SIRIUS 6 / CSI:FingerID / CANOPUS / MSNovelist / COSMIC
- Since SIRIUS 5, CSI:FingerID, CANOPUS, MSNovelist and DB search run via web service with account/licence
  (free academic, commercial via Bright Giant): https://v6.docs.sirius-ms.io/account-and-license/
- Internet off on Kaggle -> CSI:FingerID, COSMIC, MSNovelist unusable. Fragmentation trees / isotope scoring run locally
  (offline FP prediction not explicitly ruled out; unverified). Needs JVM in notebook.
- Plug-in: offline frag-tree formula score as ranker feature.

## 4. De novo / generative

All MSG de novo numbers assume the **ground-truth formula** is given (except MARLIN). We have only precursor m/z + adduct,
so feed top-k predicted formulas or formulas present in our 5 ppm pool; expect lower real performance.

### MS-GPT (arXiv 2607.23607, Jul 2026)
- https://arxiv.org/abs/2607.23607 , https://arxiv.org/html/2607.23607 . SAFE-GPT molecule LM (~100M structures pretrain),
  conditioned on MIST 4096-bit FP probabilities + formula; 100 samples (20 query groups x 5) ranked by consensus frequency; LoRA.
- MSG top-1/top-10 **23.91/28.65%**, Tanimoto 0.52/0.54, MCES 15.50/11.74; NPLIB1 29.76/41.07%. Post-audit, not audited;
  MIST checkpoint / batching exposure to the bug unverified.
- Code https://github.com/VIKI623/MS-GPT (Apache-2.0); checkpoints on SJTU Cloud
  (https://pan.sjtu.edu.cn/web/share/8f63823e4d82280636e537d1c789128d), weight licence unstated. MSG data (no NIST) — verify.
- T4: no absolute timing; est. few s/spectrum -> 400 x 100 samples ~1-2 h (unverified).
- Plug-in: de novo candidate generator into pool + ranker; consensus frequency as feature.

### FRIGID (Bohde et al., ICML 2026; arXiv 2604.16648)
- https://proceedings.mlr.press/v306/bohde26a.html , https://github.com/coleygroup/FRIGID
- 182M masked diffusion LM (GenMol-based) generating SAFE from MIST FP + formula; pretrain up to 230M structures;
  inference-time ICEBERG-guided remasking.
- MSG 18.29/22.00% (with refinement; audit-valid), base 16.09/18.19%; NPLIB1 25.03/33.37%.
- Weights Zenodo https://zenodo.org/records/19685145 (decoder, MIST encoders, ICEBERG) **CC BY-NC 4.0**; repo no LICENSE.
  Prize eligibility doubtful (forum CC BY-NC issue). Training data open (MSG, NPLIB1).
- T4: 6.58 s/spectrum base on RTX 6000 Ada -> est. 25-40 s on T4 -> 3-4.5 h for 400; refinement multiplies.
- Plug-in: candidate generator; or re-implement recipe with own weights.

### MIST + MolForge (Neo et al., NeurIPS 2025, arXiv 2508.04180)
- https://arxiv.org/abs/2508.04180 . Claimed MSG 31/40%; **audit-corrected 10.73/14.48%** (MIST bug + near-analog
  pretraining leakage). MolForge https://github.com/knu-lcbc/MolForge (licence unclear). Skip.

### DiffMS (Bohde et al., ICML 2025, arXiv 2502.09571)
- https://arxiv.org/abs/2502.09571 , https://github.com/coleygroup/DiffMS ; weights Zenodo 15122968 (MIT).
- MSG 2.30/4.25% (Tanimoto 0.28/0.39); NPLIB1 8.34/15.44%. 131 s/spectrum on RTX 6000 Ada -> infeasible on T4.

### MBGen (AAAI 2026, arXiv 2602.01643)
- https://github.com/biomed-AI/MBGen (no licence found). MSG 7.58/12.54%, NPLIB1 12.20/22.29% (as cited in MS-GPT table). Diffusion, slow.

### MS-BART (NeurIPS 2025, arXiv 2510.20615)
- https://github.com/OpenDFM/MS-BART ; weights CC BY 4.0 Figshare 30393544. MSG 1.07/1.11%; NPLIB1 7.45/10.99%. Weak.

### MADGEN (ICLR 2025, arXiv 2501.01950)
- https://github.com/HassounLab/MADGEN (MIT); data Zenodo 15036069 (CC BY 4.0). MSG 1.31/1.54%. Low priority.

### MSNovelist (Stravs et al., Nat Methods 2022)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9262714/ , https://github.com/meowcat/MSNovelist (AGPL-3.0).
  MIST+MSNovelist MSG 0/0%; NPLIB1 5.40/11.04%. Needs CSI:FingerID online. Not usable.

### Test-time-tuned LM (arXiv 2510.23746)
- https://arxiv.org/abs/2510.23746 . Spectrum+formula -> SMILES; pretrain on simulated, tune on test spectra.
  MSG top-1 3.16%, NPLIB1 12.88%. No code found. Test-time tuning idea possibly transferable to FPNet.

### FlowMS (ICLR 2026, arXiv 2603.18397)
- https://arxiv.org/abs/2603.18397 . Discrete flow matching; NPLIB1 9.15/12.05%; no MSG number seen.

### MSFlow (arXiv 2602.19912) — WITHDRAWN 2026-03-12 over possible data leakage. Ignore.

### MARLIN (arXiv 2607.04774)
- https://arxiv.org/abs/2607.04774 . Block-diffusion LM constrained by measured mass, no formula needed. Claims SOTA on NPLIB1
  among formula-free methods; no numbers in abstract, no MSG, no code. Watch.

### MS2Forge — https://github.com/huaqiyao/MS2forge (MIT), Bayesian flow networks, MSG ckpts on Google Drive; no paper/numbers (unverified).

### MS2Mol (Enveda 2023) — https://chemrxiv.org/engage/chemrxiv/article-details/6492507524989702c2b082fc ; EnvedaDark (226 NPs):
21% close match, 62% meaningful similarity. No weights. Hints host test set leans natural-product.

### Spec2Mol / MassGenie / Mass2SMILES / MolSpectLLM
- Spec2Mol MSG 0/0%, Tanimoto 0.12/0.16 (audit Table 9). MassGenie, Mass2SMILES: older, simulated-spectra-trained, no MSG numbers. Ignore.
- MolSpectLLM / SpecMol: https://arxiv.org/abs/2509.21861 , https://github.com/Eurekashen/SpecMol — Qwen2.5-7B, no MSG numbers, too heavy for T4.

## 5. CASMI context

- CASMI 2022 results: https://fiehnlab.ucdavis.edu/casmi/casmi-2022-results . 16 participants. Priority category top-1 2D:
  Nikolic_KUCA 38% (21/55), Dührkop/SIRIUS 26% (57/222); Bonus: Dührkop 29% (64/219), Nothias/Rutz 28%.
  Dührkop adduct 93-96%, formula 94-95%.
- CASMI 2026 press release (Business Wire, 2026-09-15): ~400 molecules, ~2,500 spectra (vs ~1213 we see — press figure may
  include private spectra), 2.5M training spectra / 275k molecules incl. Enveda-180; winners open-sourced.
  https://secure.businesswire.com/news/home/20260915334355/en/
- MSG leaderboard https://massspecgym.onrender.com returned HTTP 503 (unverified); MSG README keeps results under `results/`.

## 6. Open eligibility questions (verify before prize submission)
- FRIGID weights CC BY-NC 4.0; MS-GPT weight licence unstated; MolForge, MVP, MBGen, MS-BART code licences unclear.
- ICEBERG msg_all / FOAM ICEBERG: possible NIST'20 lineage through CE imputation.
- MIST bundled weights: training data unstated.

## Correction (added after cross-checking the Kaggle forum, 2026-10-17)

- FRIGID prize eligibility is **not** doubtful: the host ruled "Yes, the official FRIGID checkpoints are allowed
  for prize-eligible solutions" (Marie Killian, discussion 746430, 2026-10-13; repeated in 748882), provided the
  training code and data are referenced. Third-party mirrors (e.g. MS-GPT's bundled MIST encoders) are only
  allowed if fully reproducible (746430, 2026-10-16).
- The host also ruled that ICEBERG PubChem-atlas spectra (NIST'23-trained ICEBERG outputs) "can be used" as
  model outputs, shipped as a dataset with the generating/query code (745134).
- Spot-checked: arXiv 2606.19624 (MassSpecGym in the Wild; co-authored by David Healey of Enveda) and arXiv
  2606.29161 (GLACIER, v2 2026-10-05: MassSpecGym top-1 70.0% with contrastive finetuning, NIST'20 52.5%).
