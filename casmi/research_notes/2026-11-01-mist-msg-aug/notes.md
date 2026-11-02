# MIST msg_aug encoder probe (2026-11-01)

Offline only: nothing was submitted. Run on `pi-a100-80gb` (A100 80 GB), at about 21:30 AEDT.

## Verdict: FAIL, just below the pre-registered gate

The gate, fixed before the run: mean fingerprint Tanimoto on the 50 Enveda-180 FRIGID panel molecules must be at
least 0.40. `msg_aug` scores **0.385** (median 0.336, SE 0.022, 17 of 50 at ≥ 0.4). That is 0.015 short, inside one
SE, but the rule was set in advance, so the MIST/FRIGID line stays closed. If it is reopened anyway, that is a
post-hoc decision and should be labelled as one.

| Spectra (n) | old `mist_msg.pt` (FRIGID bundle) | `mist_msg_aug.ckpt` | `mist_msg_noaug.ckpt` |
|---|---|---|---|
| FRIGID panel, Enveda-180 (50) | 0.261 | **0.385** | 0.178 |
| check_mist Enveda-180 sample (40) | 0.247 | 0.347 | 0.179 |
| check_mist MassBank sample (40, likely in MSG train) | 0.629 | 0.729 | 0.148 |
| MSG val, official cache (19,429) | not run | 0.436 | 0.179 |
| MSG train, random 5,000 rows | not run | 0.859 | 0.191 |

The table shows means at the threshold calibrated on MSG val by MIST's own rule (matching the total on-bit count):
0.188 for `aug` and 0.110 for `noaug`. At FRIGID's 0.187, `aug` is unchanged to within 0.001. The old column was
measured at 0.187 through FRIGID's MIST code, so it compares directly.

**`noaug` is not a usable control.** It scores 0.19 on its own training spectra, against 0.86 for `aug`, read
straight from the released feature cache with no preprocessing of ours involved. The checkpoint was saved at step
30,000, against 430,965 for `aug`. So the release looks undertrained, and how much of the gain comes from the
synthetic data, rather than the new code or the longer training, is not measured.

## Setup

- **Checkpoints:** [Zenodo 23303289](https://zenodo.org/records/23303289). `mist_msg_aug.ckpt` md5
  `3e92179240e3043733627d720f4b7b86` and `mist_msg_noaug.ckpt` md5 `80594e4dfd82d15fb32a3b186e686b6b`, both matching
  the record. MSG val/train features come from `msg_compact.tar.gz` (md5 `13e86837…`).
- **Code:** [coleygroup/MIST](https://github.com/coleygroup/MIST) at `48e339be`, torch 2.14.1+cu126, Python 3.11.
  The probe script is `~/mist_probe/probe.py` on the A100, and its report is `~/mist_probe/probe_report.json`.
- **Panel:** the 50 FRIGID targets (`~/casmi-data/frigid/work/targets.json`), plus the two 40-spectrum samples
  `scripts/check_mist.py` drew (same filters, `default_rng(1)`; 40 of 40 kept in each).
- **Preprocessing:** MIST's own steps. Drop peaks above precursor + 1, sqrt-normalise, keep the top 60 above 0.001,
  assign subformulae at 20 ppm as `[M+H]+` (the old probe also used 20), keep the 10 most intense assigned peaks plus
  the CLS token, instrument `qtof`. The median is 10 annotated peaks; 2 of 130 spectra had none (CLS only, still
  scored and kept in the denominator). Truth is RDKit Morgan radius 2, 4,096 bits, no chirality.

## Not established

- **The effect on ranking.** At 0.26, FRIGID found 0 of 50 in its top 10, and top-10 Tanimoto tracked MIST Tanimoto
  (r = 0.82). Whether 0.385 gives non-zero hits was not run.
- **Use as an FPNet replacement or a candidate re-ranking feature.** Not tested.
- **T4 runtime.**
