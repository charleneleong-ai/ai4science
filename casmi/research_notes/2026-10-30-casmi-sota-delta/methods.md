# Literature delta notes, 2026-10-17 → 2026-10-30

Searched: arXiv API by submission date (MassSpecGym, tandem mass, MS/MS structure, spectrum prediction), HF
papers/models/datasets (`hf_fs search`), GitHub (`gh` commits/releases on pluskal-lab/MassSpecGym, pluskal-lab/DreaMS,
samgoldman97/ms-pred, coleygroup/*), Zenodo API, WebSearch for bioRxiv/ChemRxiv. All numbers author-reported.
arXiv dates below checked against the arXiv API on 2026-10-31; the FRIGID README table checked with `gh api`.

## New and relevant

### FRIGID / MIST leakage correction (coleygroup, 2026-10-20/21; arXiv 2604.16648 v3, 2026-10-25)
- FRIGID README: the previous MSG MIST encoder was trained with synthetic spectra of PubChem molecules chosen to
  match test-set formulae. Retrained numbers, previous struck through:
  NPLIB1 top-1 25.03% → **14.16%** (FRIGID-base 13.31%); MassSpecGym top-1 18.29% → **9.05%** (base 8.21%).
  https://github.com/coleygroup/FRIGID · https://arxiv.org/abs/2604.16648
- New MIST repo https://github.com/coleygroup/MIST (MIT, training code) and Zenodo record 23303289
  (https://zenodo.org/records/23303289, published 2026-10-20, MIT): `mist_msg_noaug.ckpt`, `mist_msg_aug.ckpt`
  (MSG + GLACIER-simulated spectra of a refined-4M molecule set; 78.45M synthetic feature rows), CANOPUS
  equivalents, `glacier_msg.ckpt`, `glacier_canopus.ckpt`, ~290 MB per MIST checkpoint. Host: allowed (751373).
- For us: explains part of FRIGID 0/50 (the bar was inflated). Our closed MIST test used the old `mist_msg.pt`
  shipped with FRIGID (experiment card: Tanimoto 0.26 on Enveda-180 vs 0.63 on MassBank). The aug checkpoint is a
  new, legal, T4-sized encoder; untested on timsTOF.

### Separability Is Not Alignability (arXiv 2610.37979, 2026-10-27, NeurIPS'26 workshop)
- https://arxiv.org/abs/2610.37979 · code https://github.com/autonlab/msms-retrieval (MIT, training code, no weights).
- MassSpecGym v1.5 retrieval, MCES split, encoders retrained to remove test-fold contamination. Frozen alignment
  (DreaMS / MS-JEPA + frozen MoLFormer) caps at MRR 0.227. Pretraining on 5.9M molecules with GLACIER + FraGNNet
  simulated spectra then fine-tuning: MRR 0.232 → 0.463; + MaxSim token scoring + candidate-aware training:
  **MRR 0.566, R@1 0.425**.
- Candidate-set negatives teach a spectrum-independent "library compound vs decoy" prior: spectrum-blind ranking
  already gives ΔMRR ≈ 0.12-0.14.
- Feasibility: training 4×A100 (paper); inference T4-fine. 3-5 days plus simulating our pools. Too big now.

### GEMS (arXiv 2610.23617, 2026-10-16)
- https://arxiv.org/abs/2610.23617. Learned bond-edit search over same-formula isomer graphs from PubChem seeds,
  guided by MIST fingerprints. MSG hit@1/10 30.2/41.3%, NPLIB1 29.6/50.8%, ~3.3 min vs FRIGID 26 h.
- No code; uses MIST and predates the leak fix; a variant reports on NIST. Not usable; watch for code.

### SPARC (arXiv 2609.35649 v2, 2026-10-18)
- https://arxiv.org/abs/2609.35649. Reference-guided test-time specialisation of a forward predictor: retrieve
  related library spectra and recalibrate fragment intensities for acquisition shift. No code found; idea only.

## Minor / no action
- MassSpecGym in the Wild v2 (2026-10-29) https://arxiv.org/abs/2606.19624; MassSpecGym leaderboard: no commits
  to main since 2026-05-08, so unchanged.
- ms-pred: one commit 2026-10-26 (HF dataset strategies), no new weights.
- MARLIN v2 (https://arxiv.org/abs/2607.04774, 2026-10-25): formula-free de novo, NPLIB1 top-1 16.9%.
- HF Hub: no new relevant checkpoints.

## Earlier 2026 papers the 2026-10-17 sweep did not list
- Selective prediction for MS/MS retrieval (https://arxiv.org/abs/2603.10950 v2): ranking-level uncertainty
  (score margin / entropy over candidates) best predicts whether rank 1 is right; fingerprint-bit uncertainty poor.
- What are we optimizing for? (https://arxiv.org/abs/2602.16507): fingerprint accuracy vs retrieval trade-off.
- MSAlign (https://arxiv.org/abs/2605.19752): frozen DreaMS + MolDeBERTa alignment (the capped route above).
