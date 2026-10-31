# Kaggle delta notes, 2026-10-17 → 2026-10-30 (pulled 2026-10-30)

Source: Kaggle MCP (`list_competition_topics`, sort New, pages 1-2: 137 topics total, 36 new since the last
sweep's 101) and `kaggle kernels list/pull`. Thread URLs are
`https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/<id>`. Forum numbers are
participant claims unless marked HOST.

## Leaderboard (CLI, 2026-10-30 ~09:00 UTC)

| rank | team | score | last sub (UTC) |
|---|---|---|---|
| 1 | Oleksandr Parfionov | 0.517 | 10-30 |
| 2 | Savko | 0.513 | 10-30 |
| 3 | Shehab Alshehabi | 0.501 | 10-30 |
| 4 | pikachu | 0.499 | 10-29 |
| 5 | CHU1PC | 0.495 | 10-29 |
| 6-7 | team051, Paweł Biernat | 0.490 | |
| 8-20 | JJ … AibePC | 0.481-0.489 | MarvinTMB still 0.484 from 10-08 |

Top moved 0.495 → 0.517 in 13 days. kunaldesale2408's notebook (§1) counted on 2026-10-27: 4,194 teams,
130 at exactly 0.472, 169 at 0.468, 190 at 0.464; gold line rank 18 at 0.477; "No member of the top ten teams
is the author of any of the 974 public notebooks"; 103 of 112 public notebooks ≥0.460 attach dmitriigluzdov's
v4g assets.

## Host rulings (new)

| thread | ruling |
|---|---|
| 750330 (HOST, sticky, 10-17) | Non-commercial open data (CC BY-NC) now allowed generally. NIST and derivatives (weights) still banned. |
| 751373 (HOST, 10-22/10-26) | Models trained only on public libraries with public code (MS2DeepScore 2.0, FIORA-OS) allowed; SIRIUS/CSI:FingerID **not** allowed; chemical language/graph models pretrained on PubChem structures allowed as priors; intermediates need generating code, not bit-for-bit (give version/download date); PubChemLite fine. **coleygroup MIST checkpoints (Zenodo 23303289, 2026-10-20, MIT, incl. GLACIER-augmented) allowed.** |
| 752297 (HOST Marie Killian, 10-26) | "Yes, automated, unsupervised processing across test spectra at inference time is allowed" — e.g. clustering test spectra, using a confidently identified test molecule to rank candidates for a similar one (confirmed again in 751373). |
| 750684 (HOST, 10-27) | Class shares will not be released; "more recent updates of the data sources used in train.parquet" allowed. |
| 752541 (HOST, 10-26) | "We reserve the right to remove submissions that are not reproducible from the leaderboard, regardless of ranking." |
| 753775, 753662, 753900 (open) | Whether participant-published checkpoints / intermediates without training code (megayak sim ranker rows, prvsiyan fp caches, ahmedberatozer pubchem-tier build, dmitriigluzdov NP knowledge table) are compliant. **Unanswered** — this is the whole public lineage, ours included. |

## Measured results posted (participant claims)

- **busyaprime 0.477 notebook submission table** (kunaldesale2408 copy, §5): leader-safe fusion 1/(3+b)+w/(3+n),
  w_top 0.25, w_mid 0.75; fragment check (Roman Tamrazov) on/off; neural model Saner Turhaner v1 vs Ahmed Berat
  Özer "long-trained" `final_model.pt`. Results: e1 0.472; long-trained neural model 0.453 (e3), 0.461 without
  fragment check (e4); NX isomer reranker weight 2: **0.414**; m1 (fragment leader moves from V2 order on top of
  0.75 order) **0.477**. → swapping in a stronger-trained neural scorer *lost* 0.011-0.019; the isomer reranker
  lost 0.058.
- **dmitriigluzdov "Data Behind" §6 versions**: removing PubChem-only candidates 0.464 → 0.409; 6 → 20 analogue
  parents 0.464 → 0.456; five fold-safe FPNets vs one: no change; his own competition-trained spectrum–structure
  re-order with leader-safe fusion: 0.459 (v27) vs 0.465.
- **Kino Lutor (750241, 10-17)**: on 130 molecules absent from candidate DBs (own split), V22 puts truth at rank 1
  only 7 times; Engine-1's analogue generator *can* produce the truth for 100/130 given suitable parents. But the
  20-parent LB test above lost 0.008.
- **blaise yonga (753400, 10-28)**: 65 biotransformation rules × 15 spectral analogues, 248 NP panel with truth
  removed: truth generated 46% / 32% / 16.5% at analogue Tanimoto caps none / 0.70 / 0.55, but class-3 MRR only
  0.080 / 0.056 / 0.030, costing 0.004-0.011 on in-pool molecules; expected net +0.004-0.010; not used.
- **Bull (751300, 10-21) and replies**: two weeks of class-3 work, local gains, no LB movement. Tra Crew:
  "a class 3 lever with consistent tiny leaderboard gains" of +0.003-0.004. Mansur: chasing recall not rerank,
  nothing clean.
- **Pathik Patel (752735, 10-26)**: failure ledger n=500 own honest split: not in any DB 21.6%, in window but
  outside top-25 14.8%, ranked 2-25 39.4%, rank 1 23.8%. Local Δ vs LB Δ: +0.0093 → +0.019; pool expansion to
  97M +0.113 → +0.003; bigger ranker +0.011 → −0.011.
- **Akhil nath reddy (753511, 10-28)**: noise model; ~125 public molecules; with 10% of molecules changing, a null
  change shows |Δ| ≥ 0.010 in 47% of comparisons; best of 7 null variants shows +0.019 public, 0 private; a
  paired offline comparison on 600 molecules has SE ≈ 0.0065.
- **Rishav Yadav (753724)**: 18M deduplicated candidate DB (COCONUT 2.0 460k + PubChem bioactive 100-1200 Da);
  no LB evidence.
- **denpugovkin Protected Tail**: adjacent same-formula swaps at ranks 2-10 when the challenger explains ≥ intensity
  in every CE context; on a 0.365 base, no score reported.
