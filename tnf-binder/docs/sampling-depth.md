# Sampling depth: what three diffusion samples could not tell us

The submission's pH-switch evidence came from **three** Boltz-2 diffusion samples per design.
This is the same measurement at **twenty** — 25 designs, 500 structures, one PROPKA run per
histidine per structure. It does not overturn the switch, but it removes two of the four
positions the campaign was built on and it changes which designs are the best ones.

Run: `run_deep.sh`, finished 2026-10-13T16:46Z. Outputs `pka_robustness_deep.csv` (62
histidine-positions) and `deep_geometry_conserved.csv` (500 structures).

## The criterion had to change, and that is the first finding

`pka_robustness.py` calls a histidine **robust** when it sits in the 6.0–7.6 window on *every*
sample. That is a maximum over samples, so it gets strictly harder as samples are added: at
n=3 one bad structure out of three is needed to fail, at n=20 one bad structure out of twenty.
Comparing "3/3" against "20/20" therefore measures the sample count, not the design.

- **7 of 62** positions are in-window on all 20 samples.
- **10 of the 20** positions that read 3/3 on three samples are not 20/20 on twenty.
- Under the all-samples rule, **0 of 25** designs carry the ≥2 switching histidines the
  submission gate requires. At n=3, most did.

The comparable quantity is the **in-window rate**, which does not depend on n. The table below
uses it, and the rate threshold is declared post hoc — it was chosen after seeing these
numbers, so the sensitivity is shown rather than one threshold quoted:

| in-window rate ≥ | designs with ≥2 such His | with ≥1 |
| --- | --- | --- |
| 1.00 (all samples) | 0 | 7 |
| 0.95 | 2 | 16 |
| **0.90** | **5** | 18 |
| 0.85 | 7 | 19 |
| 0.80 | 9 | 19 |
| 0.75 | 10 | 22 |

## Two of the four switch positions are not switches

Pooled over all 25 designs, per position:

| His | n | mean pKa | median spread | max spread | mean in-window rate |
| --- | --- | --- | --- | --- | --- |
| **67** | 13 | 6.82 | 1.56 | 2.45 | **0.94** |
| 72 | 11 | 6.30 | 1.98 | 3.37 | 0.84 |
| 70 | 12 | 6.05 | 0.45 | 1.39 | 0.83 |
| 38 | 7 | 5.82 | 1.89 | 3.14 | 0.56 |
| **21** | 13 | 6.20 | **4.28** | 5.73 | 0.56 |
| **35** | 6 | 5.90 | **5.70** | 6.73 | 0.42 |

**His21 and His35 have no determined pKa.** A median spread of 4.3 and 5.7 pH units across
structures the folding model scores as equally good is not a measurement with a window in it —
the predicted pKa is a property of the sample, not of the design. Three samples hid this: on
three structures His21 in `helix0_R21H_E70H` read a spread of 0.50 and was called robust; on
twenty it spans 5.21.

His67 is the one position that holds up (rate 0.94, spread 1.56), with His70 tight but centred
low (6.05, i.e. near the bottom of the window) and His72 intermediate.

Every design that passes ≥2 positions at rate ≥0.90 does so on **E67H plus E70H or R72H**.
Where R35H is also present it contributes nothing to the pass.

## Geometry over the same 20 samples

No design is clash-free on a majority of its samples; the best is 9 of 20.

| design | His ≥0.90 | clash-free | min area Å² | min focus | median area Å² |
| --- | --- | --- | --- | --- | --- |
| **helix0_R21H_E67H_E70H** | 2 | **9/20** | 580 | 0.298 | 746 |
| helix0_E67H_R72H | 2 | 7/20 | 553 | 0.224 | 797 |
| helix0_R35H_E67H_E70H *(submitted)* | 2 | 7/20 | 552 | 0.165 | 714 |
| helix0_R35H_E67H_R72H | 2 | 5/20 | 455 | 0.254 | 723 |
| helix0_E67H_E70H | 2 | 3/20 | 431 | 0.090 | 786 |

All 25 designs, both axes, are in
[`deep_geometry_conserved.csv`](deep_geometry_conserved.csv) and
[`pka_robustness_deep.csv`](pka_robustness_deep.csv).

`epitope_focus` here is the share of buried target area on the conserved TNFR2 site, which is
what `submit.py` declares. It is not comparable with a focus computed against the design pocket
(113/115/143/144) — on one submitted design the same structure reads 0.09, 0.18 and 0.56 under
the three epitope definitions the campaign has used. See
[`scripts/score_designs_geometry.py`](../scripts/score_designs_geometry.py).

## What this does and does not establish

- It **does** establish that R21H and R35H cannot be claimed as pH switches, and that the
  ≥2-robust-histidine gate as written is unsatisfiable at this sampling depth.
- It **does not** refute the switch at E67H/E70H/R72H. Those rates are high and their spreads
  are small enough that the window means something.
- It is still PROPKA on predicted structures. No experimental pKa, no bound-state measurement,
  and no pH-dependent binding assay is involved anywhere in this.
- The 0.90 rate threshold is post hoc. The pre-registered criterion was "every sample", which
  at n=20 selects nothing.
