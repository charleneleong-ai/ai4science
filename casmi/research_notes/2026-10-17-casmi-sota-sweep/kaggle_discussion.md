# Kaggle discussion and public-notebook notes (pulled 2026-10-17)

Source: Kaggle internal API (`discussions.DiscussionsService/GetTopicListByForumId`, forum 11293906; 101
topics, all fetched with comments) and `api/v1/kernels/pull`. Topic index by votes:
[`kaggle_topic_index.txt`](kaggle_topic_index.txt). Thread URLs are
`https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/<id>`.
Forum posts are participant claims unless marked HOST/ADMIN; numbers are self-reported and unverified.

## Competition facts (API metadata + host posts)

- 400-ish hidden molecules, `leaderboardPercentage: 31` → public LB is ~31% (~124-132 molecules; one molecule
  from absent to rank 1 ≈ +0.008). Deadline 2026-12-14, merger deadline 2026-12-07, 5 subs/day, 540 min GPU.
  3,463 teams.
- LB top 2026-10-17: CHU1PC 0.495, pikachu 0.495, Oleksandr Parfionov 0.493, MarvinTMB 0.484 (set 2026-10-08,
  before the public 0.44 lineage existed), Oliver Ford 0.468, then a dense band 0.456-0.467.
- Host, 741359: "We are especially excited to be providing spectra for many molecules whose structures have not
  been previously described."
- Host, 745758 (Marie Killian): will not give the class breakdown; "our hope is that there are enough novel
  structures to make it worth working on them."
- Host, 742274 (David Healey): class 3 is absent from **both PubChem and COCONUT**. Empty guesses are dropped.
  Metric: https://www.kaggle.com/code/metric/casmi-mean-reciprocal-rank
- Host, 748882: all three classes in both public and private; live scorer = published scorer v13; appending
  distinct valid guesses can never lower a row. Electron-mass correction applied **only to the hidden test**,
  not to Enveda-180 in train.parquet nor the Zenodo release: "subtract 0.000548 Da from positive-mode m/z values"
  to align them (see also 743395, 746954).
- Admin, 745715: no ID overlap between train/test or dummy/rerun test; "relatively low-observation test set …
  careful that they don't interpret noise with causation."
- Organisers (press release via Yahoo Finance mirror of BusinessWire, 2026-09-15): "The molecules are real or
  potential natural molecules and their analogs"; ~400 molecules, ~2,500 spectra, unpublished, Enveda
  isolation/profiling platform; CASMI founders Schymanski and Neumann support.
  https://finance.yahoo.com/healthcare/articles/enveda-launching-next-alphafold-contest-120700676.html
- Hidden test: timsTOF only, 157-1,159 Da neutral mass, 1-16 spectra per molecule (data page, as quoted in
  dmitriigluzdov notebook).

## Eligibility rulings that open new routes (host answers)

| Resource | Ruling | Thread |
|---|---|---|
| ICEBERG PubChem atlas (91.4M PubChem structures, predicted by NIST'23-trained ICEBERG) | "The first reading applies. The atlas spectra are model outputs publicly offered by the authors and can be used." Shipping only the predicted spectra for your candidates as a dataset, with the generating/query code, is fine | 745134 |
| FRIGID official checkpoints (Zenodo 19685145, CC BY-NC) | "Yes, the official FRIGID checkpoints are allowed for prize-eligible solutions" with reference to training code/data | 746430, 748882 |
| TabPFN v2 / 2.5+ | eligible | 748882 |
| PubChem CID-PMID / CID-Patent / CID-SID counts | allowed as ranking features | 748882 |
| LOTUS, NPAtlas | acceptable structure sources | 748882 |
| DreaMS SSL backbone; GLACIER/ICEBERG (MassSpecGym); MassSpecGym-trained models | allowed | 748270 (table), 745604, 744338 |
| CFM-ID 4; GNPS-NIST14-MATCHES; newer GNPS/MassBank/MoNA; Enveda-180 Zenodo (816 extra structures) | allowed | 748270, 745148 |
| NIST and anything trained on it | disqualified | 744470 |
| Third-party fine-tuned weights without public training code | **not allowed** | 748882 |

## Class-share estimates (participant LB probes; disagree)

- prvsiyan analog-propagation notebook via haideptry (741597): ~16% class 1 / 45% class 2 / 39% class 3,
  derived from a library-only submission scoring 0.151 with class-1 MRR ≈ 0.93.
- Bull (748524): fitted from several submissions, ~28% / 20% / 52%.
- Antriksh Arya (748029): "about 132 public molecules and almost half of it are novel".
- haideptry (741745): tier 1 10-15%, tier 2 45-55%, tier 3 30-40%.
- MarvinTMB (746302, 0.484 at the time): analytical ceiling ~0.76, "realistic maximum probably around 0.55";
  "You do not need NIST or any additional data for a better score" (744470).

## What moved / did not move the public LB (single submissions; noise ±0.008/molecule, ±0.006 seed)

From dmitriigluzdov "From Spectra to Structures" §13-14 (133 votes, run 2026-10-17):
https://www.kaggle.com/code/dmitriigluzdov/casmi-26-from-spectra-to-structures
- v49 0.446 (rel>150) → v51 **0.450** deeper trees (depth 7) in second-engine rankers (from
  ghazarosbarseghyan91) → v52 **0.451** standardised ChEBI/LIPID MAPS (sergiomalv) in engine 2 → v53 more
  seeds (8+6), score not shown. So the public best is now **0.451**, not 0.446.
- Did not help on the v43 (0.440) pipeline: GLACIER per stepped-CE step merged FIORA-style (0.436); MSEntropy
  centroiding before GLACIER (0.440); separate merged spectra per polarity/adduct in engine-2 fingerprint branch
  (0.424); **RRF with JESTR spectrum-structure similarity (0.418)** — released JESTR encoders are NPLIB1-trained;
  **CFM-ID 4 spectra of the LOTUS in-silico DB**, same-formula reorder for library sim < 0.7 (0.440 at w=0.25,
  0.422 at w=1.0); fragment score joined to ICEBERG where GLACIER missing (0.441 → 0.440).
- Aryan Awasthi (749908): what separates 0.446+ from 0.443 is rel>150 promotion, metric-exact engine-2 rank-2
  insert and the CC-BY ChEBI/LMSD set; unpinned Docker image cost 0.014 (0.441 → 0.427).
- Antriksh Arya (749412): masked probes — on E006 disagreements, E006 #1 right ~13, engine-2 #1 right ~3;
  engine-2 top-1 inserted at rank 2: +0.002 on two bases.
- Aryan Awasthi (749251): an independent spectrum→fingerprint MLP + ranker over PubChem/COCONUT merged by RRF,
  top-1 protected: 0.437 → 0.438. On his class-2-like split (490 enveda-180 structures without public spectra),
  fragmentation features in a LambdaRank ranker 0.502 → 0.627, analog features → 0.660, but LB moved only
  +0.03-0.04 for those steps. Popularity puts the answer #1 for 67% of enveda-np-examples but ~0% of enveda-180
  no-public-spectra molecules (748068).
- Udam Liyanage (743254): ~98% of database-route errors are same-formula isomer errors; offline gains < 0.03
  invisible on LB; generated-candidate insertion "has not yet paid off on the leaderboard, even when only
  inserted for molecules predicted to be novel"; deleting the truth from the DB overstates generator reach ~3×.
- NeckBeard (743254): MetFrag-lite tweaks held out: +0.023 ± 0.007 panel MRR (2-bond cleavage ×0.6, linear
  intensity, 0.005 Da, H shift −2..+3); LB 0.341 → 0.336 (nothing).
- Bull (748524): ICEBERG re-rank +0.02 on class 2 overall but ~0 on structures not in MassSpecGym (training
  exposure); fixed-rank PubChem insertion +0.01-0.02 locally, LB 0.295 → 0.288.
- Victor alexandre (742088): perfect-rerank ceiling of the old 0.324 pipeline's top-25 ≈ 0.40; filler/junk
  SMILES zeroed a submission.
- DancingLumberjack (742088): PubChem as full candidate source +0.002 on LB.
- haideptry (741745): blind PubChem isomer expansion hurt; ±8.5 ppm better than wider.

## Hints about leaders' methods

None of the top-4 teams has posted a method. Indirect signals only:
- MarvinTMB reached 0.484 on 2026-10-08 independently of the public E006 lineage and says no extra data is
  needed (unverified).
- hengck23 (745070, 744343, 744893): peak → fragment-formula annotation → structure ("hard SIRIUS to neural
  SIRIUS"); train a forward SMILES→MS/MS model ("old experienced kaggler will know what this is for" — i.e.
  synthetic training data / forward-model rescoring); DreaMS kNN as a coverage probe.
- Adarsh (745468): formula correct in >95% of validation, failures are isomers; next step fragment-formula
  plausibility.
- Host tutorial (inversion/casmi-denovo-tutorial-notebook, 358 votes): a 6+6-layer spectrum→SMILES BPE
  transformer trained on non-Enveda-180 train, 25 samples per spectrum — the host's own pointer to de novo.

## Other public resources worth knowing

- flexonafft "Closing the Evidence Gap" (2026-10-17) — union of old prvsiyan BIO set and sergiomalv BIO set in
  engine 2 on the "V22 0.451 recipe"; score unknown.
- hengck23 DreaMS dataset (already attached in our fork); Aryan Awasthi validation sets:
  https://www.kaggle.com/code/aryanawasthiteykik/casmi26-validation-that-behaves-like-the-test and
  https://www.kaggle.com/code/aryanawasthiteykik/casmi26-gpu-fingerprint-fragmentation-ranker
- dmitriigluzdov identity-disjoint folds dataset (identity and scaffold folds).
