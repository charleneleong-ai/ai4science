# AI4Science learning log

Dated, distilled entries — newest first. Each entry: what it is / why it matters / how it transfers to my work / what to learn next. Maintained by the `ai4science-log` Claude skill (or by hand, same format). Evergreen overview + reading list live in [`README.md`](README.md#for-my-own-learnings--an-intro-to-ai4science-modelling).

<!-- NEW ENTRIES GO DIRECTLY BELOW THIS LINE (newest first) -->

## 2026-10-07 — CASMI 2026: no public database covers the test chemistry, because it isn't natural products

**Source:** measured — [COCONUT 2.0](https://coconut.naturalproducts.net/download) (CC0, 474k structures after collapsing to `inchikey14`) against the competition data; writeup in [`casmi/docs/2026-10-07-casmi-spike.md`](casmi/docs/2026-10-07-casmi-spike.md) (grounded 2026-10-07)

**What it is:** COCONUT covers **0 of the 400** public test molecules, and 8% of held-out train molecules. The test answers were recovered by exploiting the leak — a hash join on `(adduct, precursor, n_peaks)` resolved 400/400, independently confirming the duplication found earlier by cosine. All 400 come from one library, `enveda-180`; 36% contain a halogen; the structures are brominated cyclopropyls, triazoles, spiro-morpholines.

**Why it matters:** despite the CASMI name and Enveda being a natural-products company, the test set is **synthetic drug-like chemistry**. So a submission needs PubChem-scale retrieval, not a natural-product database, and that is the dominant unsolved cost of the whole task. It also kills the reason I had given for ICEBERG — see the correction in the entry below.

**How it transfers:** the generalisable habit, and the third time today it paid: *check what the data is before taking the domain label at face value*. "CASMI" plus "Enveda" implied natural products, the literature reinforced it, and a reasonable inference was wrong. Same shape as [[touchstone]] today, where the spec asserted an N₂S₂ motif the verifier was never scoring. A label is a hypothesis about data, and it is cheap to test.

**To learn next:** whether the private split is drawn like the public one (presumably, but unverifiable); and what a PubChem-scale pool does to a candidate window currently at ~46 — that number decides whether the 0.70 ranking result means anything outside train's own chemical space.

## 2026-10-06 — CASMI 2026: the public test set is duplicated in train, so the leaderboard measures lookup

**Source:** measured directly on the competition data — [`train.parquet`](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/data) (2,539,608 spectra / 275,810 unique `inchikey14`) vs `test.parquet` (1,213 spectra / 400 molecules), 2026-10-06

**What it is:** Every one of the 400 public test molecules has a train spectrum at peak-cosine ≥0.95; the median best cosine is exactly **1.000**. Spot-checking five molecules against the full train scan: cosine 1.0000, identical peak counts (183/183, 219/219, 244/244, 31/31, 287/287), identical m/z values, all from the `enveda-180` library — and each match carries `normalized_smiles` and `inchikey14`, i.e. the answer. The public test spectra are duplicates of training rows.

**Why it matters:** the public leaderboard is an **exact-lookup task**, not structure elucidation. Join test→train on precursor + peaks and read off the structure. This explains the ~0.87 MRR@25 figures in public repos — and inverts their meaning: they are evidence the join works, not that a model works. It also makes the public LB useless for model selection, because every comparison on it measures lookup ability. Only a molecule-disjoint split of train predicts private-set performance, which is why `inchikey14` already being a column matters: it *is* the metric's equivalence class, so the split is a groupby rather than a chemistry problem.

**How it transfers:** this is the dataset-level version of the thing [[touchstone]] kept teaching today — *the measurement is only as good as the object it is measuring*. There, a 0.02 Å parse cutoff silently scored N₂S₁ while the spec claimed N₂S₂; here, a leaked duplicate silently scores lookup while the leaderboard claims elucidation. Both are verifier-integrity failures, not model failures, and both were invisible until something was actually run. The habit that catches them is the same: before trusting a number, establish what it is a number *of* — negative controls (random-pair cosine mean 0.0203, max 0.107) and self-controls (1.0000) before believing a headline.

**The trap I nearly fell into:** my first verification sampled only early train batches, found 40/40 precursor-window matches were *different* molecules, and appeared to refute the leakage. Both results were true — most same-mass matches are different structures, while one identical spectrum exists somewhere in 2.5M rows. A verifier that examines the wrong subset produces a confident wrong answer, so scope the control to the same search space as the claim.

**To learn next:** whether the private split is also leaked (if so the competition measures nothing, which seems unlikely); and the real dependency this exposes — a train-only candidate pool makes held-out structures absent by construction, so a realistic retrieval pool needs an external structure DB (PubChem/COCONUT). Until then the harness can measure ranking-given-truth-in-pool, but not recall.

## 2026-10-06 — MS/MS structure elucidation: the generator moved, the verifier is still the bottleneck

**Source:** [MS-GPT 2607.23607](https://arxiv.org/abs/2607.23607) · [MARLIN 2607.04774](https://arxiv.org/abs/2607.04774) · [FlowMS 2603.18397](https://arxiv.org/html/2603.18397v1) · [DiffMS 2502.09571](https://arxiv.org/html/2502.09571v2) · [de-novo review, Molecules 31(5):769](https://www.mdpi.com/1420-3049/31/5/769) · [ms-pred/ICEBERG](https://github.com/coleygroup/ms-pred) (grounded 2026-10-06, for [Enveda CASMI 2026](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra))

**What it is:** Predicting a molecule's 2-D structure from its tandem mass spectrum. Two families: *retrieval* — predict a fingerprint from the spectrum and rank database candidates (CSI:FingerID/SIRIUS) — and *de novo* — generate the structure outright. De novo has moved fast through three eras: fingerprint-conditioned RNNs (MSNovelist) → end-to-end sequence models → **graph-native diffusion under formula constraints** (DiffMS, FlowMS), and now molecule-language posterior querying (MS-GPT: 29.8%/41.1% top-1/top-10 on NPLIB1, 23.9%/28.7% on MassSpecGym). MARLIN drops the ground-truth-formula crutch entirely.

**Why it matters:** I had assumed retrieval+rerank was the winning shape and de novo a niche for novel structures. Half wrong. De novo is where the research frontier moved — graph-native models were the first to break 0% exact-match on leakage-controlled MassSpecGym. But the *competition*-relevant claim survives and sharpens: **candidate recall is ~solved (≈100% within a mass window) and ranking is the whole problem.** The published winning CASMI shape is a three-class hybrid — library search, fingerprint retrieval, de novo for the tail — not any one family. So "generator vs verifier" doesn't map onto "de novo vs retrieval"; the generator is candidate *supply*, and the verifier is the ranker, whichever family supplied them.

**How it transfers:** This is the verification-first thesis in a domain with a free, cheap, *physics-adjacent* verifier — a forward model. ICEBERG predicts a spectrum from a candidate structure, so you can score any candidate by round-tripping it back to the observed spectrum, exactly the recursive check in the [[touchstone]] stack (generator proposes a site → independent oracle judges it). ICEBERG 2.1 (2026-07-07) ships a GPU-fast pretrained NIST'23 model, so it's a reranker you rent rather than train. **Correction (2026-10-07):** I also argued for it because it beats MassFormer *specifically on natural products* (0.627 vs 0.568 cosine) "which is CASMI's domain". Measured on the data, that premise is wrong — all 400 public test molecules come from `enveda-180`, 36% carry a halogen, and COCONUT covers 0 of them. The test set is synthetic drug-like chemistry, so the natural-product advantage is irrelevant here and NIST'23's licensing needs checking before any plan leans on it. The OOD-switching instinct also ports: when no library neighbour exists, fall back from cosine to fragment-explainability and de-novo candidates — a when-to-trust-imagination crossover, same as the [[touchstone]] cutoff lesson that a 0.1 Å parse choice silently changes which object you are scoring.

**The number that matters, and a caution:** public-leaderboard scores are reportedly ~**0.33** MRR@25, against **0.87** self-reported in public repos on 50-query subsets. Treat every repo README figure as unanchored until reproduced — the same trap as the BoltzGen Cu baseline whose provenance evaporated. *(The 0.33 is from a search summary of a third-party repo, not the leaderboard; unverified.)*

**To learn next:** reconcile a real contradiction in the sources — the review reports SOTA at ~4.1% top-10 on leakage-controlled MassSpecGym while MS-GPT claims 28.7%; one of those is a different split or metric, and knowing which decides whether de novo is worth a CASMI arm at all. Then: does an ICEBERG rerank actually move MRR@25 over spectral cosine, or is the forward model's error larger than the ranking signal?

## 2026-06-16 — Boltz-2 (co-folding + binding affinity)

**Source:** [paper (PubMed 40667369)](https://pubmed.ncbi.nlm.nih.gov/40667369/) · [boltz.bio/boltz2](https://boltz.bio/boltz2) · reliability eval [arXiv 2603.05532](https://arxiv.org/abs/2603.05532) (grounded 2026-06-16)
**What it is:** First open biomolecular co-folding model (MIT + Recursion, Jun 2025) to *jointly* predict 3D complex structure **and** binding affinity — near-FEP accuracy on the FEP+ benchmark at ~1000× lower cost, trained on ~5M affinity measurements; beat all CASP16 affinity entrants.
**Why it matters:** It's the cheap in-silico *verifier* I proposed for the MaterialHack protein–metal stack — an independent judge of a designed binder. But the non-obvious catch: a 2026 reliability study found its affinities aren't a trustworthy *quantitative* predictor vs. rigorous physics (ESMACS) on some compound sets, and InteractBind ([2605.24045](https://hf.co/papers/2605.24045)) shows co-folding models nail binding *likelihood* yet miss binding-*site* localisation.
**How it transfers:** The recursive verification point — *even the verifier needs verifying.* Use Boltz-2 confidence (ipTM/PAE) + held-out real data as the judge, but treat its affinity as a *ranking* signal, not ground truth; calibrate / flag OOD rather than trust the number. Generator = my binder designer; verifier = Boltz-2 co-fold; meta-verifier = held-out assay data + a physics spot-check.
**To learn next:** how Boltz-2's confidence heads behave OOD (metal-coordinated systems aren't its training sweet spot); whether to ensemble it with a physics binding-energy estimate.

## 2026-06-16 — Orb (universal interatomic potential)

**Source:** [Orb v1 — arXiv 2410.22570](https://arxiv.org/abs/2410.22570) · [Orb-v3 — arXiv 2504.06231](https://arxiv.org/abs/2504.06231) (grounded 2026-06-16)
**What it is:** A foundation model for atoms — input a 3D structure, output energy + forces; a fast drop-in for DFT in geometry optimisation, MD, Monte Carlo. Recipe: diffusion/denoising pretraining → supervised NNP. v3 charts an equivariance/conservatism/sparsity Pareto frontier and reaches the mesoscale (>10k atoms).
**Why it matters:** It's the "fast learned surrogate for quantum reality" — the verifier/world-model layer of the materials stack, and Orbital's crown jewel.
**How it transfers:** Same surrogate-in-a-loop pattern as my world-model and forecasting work. My entry point is the trust layer (UQ/OOD/conformal) that decides when to trust Orb vs. recompute DFT. Learn: conservative vs non-conservative forces, why energy conservation matters for MD stability.
**To learn next:** equivariant architectures (MACE/e3nn); how `equigrad` works.

## 2026-06-15 — When to Trust Imagination (World-Action Models)

**Source:** [arXiv 2605.06222](https://hf.co/papers/2605.06222) (grounded 2026-06-15)
**What it is:** World-Action Models imagine future frames + actions but execute a fixed chunk blindly. Reframes execution as a *future–reality verification problem*; adds FFDC, a lightweight verifier that compares predicted vs. real observation each step and replans when they diverge.
**Why it matters:** They didn't improve the world model — they bolted a *verifier* on the same imagination, and real-world success went 45% → 80%. Cleanest single instance of "the verifier is the asset."
**How it transfers:** Direct analog of RLVR (reward = verifier) and conformal forecasting (band = verifier). The materials version is a trust layer on Orb; the harder twist there is that ground truth (DFT) is expensive, unlike a robot's free camera frame.
**To learn next:** how FFDC trains its verifier on synthetic negatives — does that idea port to UQ for surrogates?
