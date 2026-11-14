# Kaggle forum and notebooks: scorer evidence since 2026-10-30 (pulled 2026-11-14)

Delta on [2026-10-30 notes](../2026-10-30-casmi-sota-delta/kaggle_discussion.md) and
[2026-11-01 leaders](../2026-11-01-leaders/notes.md). Pulled through the Kaggle MCP (topics list, topic
reads, content/dataset search, notebook source read as text, never executed). All LB numbers are single public
submissions reported by participants (public SE about 0.014; one notebook reports ~0.002 run-to-run drift from
time guards). None verified by us.

## Field state

- Leaderboard: MCP `get_competition_leaderboard` returned Unauthenticated, so the top is taken from the brief
  (leaders 0.540). Not verified here.
- Best public notebooks now show **0.488** (alexchilton
  [casmi26-lb-0-488-peaks-to-structures-network](https://www.kaggle.com/code/alexchilton/casmi26-lb-0-488-peaks-to-structures-network),
  ozertuu [casmi26-from-peaks-to-structures](https://www.kaggle.com/code/ozertuu/casmi26-from-peaks-to-structures),
  heliosli [casmi26-0-488](https://www.kaggle.com/code/heliosli/casmi26-0-488)). Dozens sit at 0.487. All are the same
  lineage as our fork (public scorer + popwenman v80 z-blend, extlib rank-1 overlay).
- Forum: 159 topics total; about 25 new since 10-30. None from a top-10 team describes a method.

## Host rulings since 10-30 relevant to a new scorer

| Ruling | Thread | Consequence for scorer options |
|---|---|---|
| "Officially released DreaMS checkpoints are fine to use" (David Healey, 11-05). Follow-up asking whether this covers `ssl_model.ckpt` (Zenodo 10997887, CC BY 4.0) as a frozen encoder + head trained on comp data is **unanswered** | [754610](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754610) | DreaMS fine-tune is allowed |
| "Must a final submission be reproducible end-to-end from the competition data and openly documented resources? Yes." | 754610 | Any new scorer needs its training code published or packaged |
| "MassSpecGym and MIST models are fine" (11-11), answering whether MSG (NIST used only for cleaning) and MIST `msg_noaug` are allowed. The GNPS-NIST14-MATCHES follow-up is unanswered | [757307](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/757307) | **MassSpecGym is cleared as training data**, so MSG-trained checkpoints (JESTR, MVP, FLARE, GLACIER/ICEBERG-MSG, MIST) are usable if their data is MSG-only |
| EnvedaDark: "IIRC it was not in fact published"; would follow normal external-data rules | [754732](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/754732) | No timsTOF NP training set from Enveda |
| `collision_energy_ev` is the nominal value; a list means multiple acquisitions merged into one spectrum (11-05) | [753419](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753419) | CE-conditioned forward models (FIORA, GLACIER) get nominal eV, not the PASEF ramp |
| Earlier, still binding: MS2DeepScore 2.0, FIORA-OS, MIST Zenodo 23303289 allowed; SIRIUS/CSI:FingerID not; CC BY-NC allowed; NIST-derived banned | [751373](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/751373), [750330](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/750330) | |

Open eligibility questions: [756893](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/756893)
(shipped ranker tables without generator), [758233](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/758233)
(small index lists), [757309](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/757309)
(PubChem/COCONUT index redistribution).

## Measured scorer swaps on the public LB (participant-reported)

Main source: the results tables reproduced in zz3t3rn41/guhongbin-lineage notebooks, read from
youyouwou's [casmi26-v225v589np](https://www.kaggle.com/code/youyouwou/casmi26-v225v589np) source (the header table is
the account owner's own submissions, "one notebook per row", with submission ids for some rows), and the bgunshui
0.480 notebook text embedded in ozertuu's
[0-500-casmi26-dual-neural-physical-isomer](https://www.kaggle.com/code/ozertuu/0-500-casmi26-dual-neural-physical-isomer).

| Change | Base | Result | Δ |
|---|---|---|---|
| bgunshui: 0.5·z(public) + 0.5·z(self-trained kit retrain, `bgunshui/casmi26-scorer-train-public`, ~2.3 h T4) on ranks 2-25 | 0.479 | 0.480 | +0.001 |
| self-trained scorer **replacing** the public one | 0.479 | 0.460 | −0.019 |
| bgunshui ensemble with popwenman v80 as the second model | 0.480 | 0.480 | 0 |
| second-scorer weight 0.5 → 0.3, 0.4 or 0.6 | 0.480 | 0.474 each | −0.006 (identical for all three weights, which is odd; possibly a pipeline fallback) |
| learned isomer reranker (`V8` RR or BOTH) on Gold Apex V22 | 0.481 | 0.460 / 0.461 | −0.020 |
| romanrozen v9 `V8='BOTH'` + bgunshui scorer v2 | 0.481 (page) | 0.461 | −0.020 |
| retrained FPNet / fe_D, 7 separate swaps (offline +26.5% Tanimoto vs public) | 0.456 | 0.418-0.443 | all lost |
| attaching popwenman v80 scorer kernel (+ dropping LOTUS backfill, + docker image; not separated) | 0.484 | 0.487 | +0.003 (confounded) |
| youyouwou v211 "ftchallenge": NP-fine-tuned v8-family scorer (`youyouwou/casmi26-v9np-bundle`) may promote one same-formula candidate to rank 1 when its z-margin ≥ 0.8 and library match < 0.90 | 0.487 | **0.424** | −0.063 (the run also reorders scorer discovery to prefer `scorer-v2-model`, so not a clean single change) |
| youyouwou v213 "ftbase": test-set molecular-network rerank (λ 0.25, τ 0.6) | 0.487 | 0.479 | −0.008 |
| youyouwou v220 "triple": `POP_W` 0.12 → 0.25 (not a scorer change despite the name) | 0.487 | 0.477 | −0.010 |
| FIORA-OS negative-mode forward-model blend, weight 0.25 (offline panel MRR only) | 0.9583 | 0.9444 | −0.014 offline |
| test-set network propagation, neighbour weight 0.25 (offline panel) | 0.9706 | 0.9674 | −0.003 offline |
| our own (2026-11-06): public + 3 seed retrains, equal z | 0.479 | 0.463 | −0.016 ([notes](../2026-11-03-seed-retrains/notes.md)) |

Reading: every scorer change in this lineage is within ±0.003 when it is a 50/50 z-blend, and loses 0.01-0.06 when the
new model gets more say (replace, re-weight, rank-1 challenger). Nobody has posted an LB number for a scorer of a
*different family* (MIST, DreaMS fine-tune, JESTR/MVP/FLARE, MS2DeepScore vs library, FIORA cosine) blended into
the 0.487 base.

## Public scorer retrains now on Kaggle (all kit-architecture unless stated)

All derived from `sanpier/casmi-spectrum-model-v1` (Apache-2.0 kit), so they are near-clones of what we already blend.

| Dataset | What it says | Licence |
|---|---|---|
| [bgunshui/casmi26-scorer-v2/v4/v5/v7/v8-model](https://www.kaggle.com/datasets/bgunshui/casmi26-scorer-v8-model) | successive retrains, v2 used by 59 kernels | CC0 / unknown |
| [godgod3/casmi26-scorer-train-v80-1](https://www.kaggle.com/datasets/godgod3/casmi26-scorer-train-v80-1) | copy of an 80-epoch run | unknown |
| [youyouwou/casmi26-v9np-bundle](https://www.kaggle.com/datasets/youyouwou/casmi26-v9np-bundle), v9large, v9b | v8-family fine-tuned on natural-product libraries (per the v211 code comment); no description | CC0 |
| [harshitkmr/casmi-ours-scorers-b](https://www.kaggle.com/datasets/harshitkmr/casmi-ours-scorers-b) | b160w (width 288, 4 layers, best epoch 99, kit diagnostic MRR@25 0.725) and b80s; train.parquet only; training script included | Apache-2.0 |
| [sardorishiguro/casmi26-scorer-v1130-model](https://www.kaggle.com/datasets/sardorishiguro/casmi26-scorer-v1130-model) | title claims "Ours 0.6603 MRR" (metric/split unstated) | CC0 |
| [studyexchange/casmi-spectrum-model-ours-s26/s27](https://www.kaggle.com/datasets/studyexchange/casmi-spectrum-model-ours-s27) | kit retrained with the real per-spectrum adduct index (train/serve agreement fix) | CC0 |
| [czhuangsjtu/casmi26-scorer-v2-phys](https://www.kaggle.com/datasets/czhuangsjtu/casmi26-scorer-v2-phys) | "Physics-Informed Dual Encoder", 127 MB, no description | CC0 |
| [prvsiyan/casmi26-self-trained-scorer-seeds](https://www.kaggle.com/datasets/prvsiyan/casmi26-self-trained-scorer-seeds) | seed retrains | CC BY-NC-SA |

No LB number found for harshitkmr b160w, studyexchange real-adduct, czhuangsjtu phys, or sardorishiguro v1130.

## Other new threads (not scorer-specific)

- [756504](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/756504) Yohann Roux,
  open pipeline (0.313, [GitHub](https://github.com/Majyles/kaggle-casmi-2026)): class 1 is cross-instrument (median
  entropy similarity ~0.90 timsTOF vs other instruments, against ~0.99 between public libraries); simulate class 1 by
  hiding the whole instrument family. Errors are "almost all same-formula isomers". COCONUT median nearest-neighbour
  Tanimoto ~0.81, so class-3 validation must remove neighbours too.
- [755319](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/755319): fork
  saturation; dauriel notes forks "rise and fall together on the private board".
- Jake (decotoj, in 756504) reports a ground-up non-public-lineage pipeline at 0.358.
- [753724](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra/discussion/753724): 18M-structure
  candidate DB (COCONUT + PubChem bioactive), not a scorer.
