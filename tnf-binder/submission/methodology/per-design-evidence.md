# Per-design evidence

Ranked on **ΔQ floor** — the charge the interface gains going pH 7.4 → 6.0, with every
histidine scored on whichever of its three predicted structures gives the *smaller*
swing. A lower bound, not a best case. One ideally-placed histidine gives at most ~0.62.

`his_robust` is how many of the design's histidines stay inside 6.0–7.6 on **all three**
structures. 8 of 23 across the set do not — those designs are demoted, not dropped.

| # | design | aa | His | ΔQ floor | robust | pKa (design structure) | buried Å² | cross-react | ipTM | Cys | novelty |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tnf2_01 | `helix0_R21H_E67H_E70H` | 91 | 3 | **1.38** | 3/3 | 7.05 7.05 6.43 | 748 | 0.64 | 0.153 | 1 | 3/4 |
| tnf2_02 | `helix0_R21H_E67H_R72H` | 91 | 3 | **1.14** | 2/3 | 7.10 7.12 6.65 | 777 | 0.8 | 0.162 | 1 | 3/4 |
| tnf2_03 | `helix0_E67H_R72H` | 91 | 2 | **1.03** | 1/2 | 7.12 6.66 | 528 | 0.7 | 0.159 | 1 | 3/4 |
| tnf2_04 | `helix0_R21H_E70H` | 91 | 2 | **0.99** | 2/2 | 7.04 6.44 | 889 | 0.7 | 0.175 | 1 | 3/4 |
| tnf2_05 | `helix0_E67H_E70H` | 91 | 2 | **0.89** | 2/2 | 7.05 6.44 | 716 | 0.74 | 0.151 | 1 | 3/4 |
| tnf2_06 | `free1_T73H` | 74 | 1 | **0.37** | 0/1 | 6.27 | 923 | 0.79 | 0.634 | 0 | 4/4 |
| tnf2_07 | `free1_parent` | 74 | 0 | **0.0** | n/a | — | 802 | 0.77 | 0.497 | 0 | 3/4 |
| tnf2_08 | `helix0_R35H_R72H` | 91 | 2 | **0.73** | 1/2 | 6.07 6.66 | 694 | 0.89 | 0.216 | 1 | 3/4 |
| tnf2_09 | `helix0_R72H` | 91 | 1 | **0.57** | 1/1 | 6.66 | 600 | 0.78 | 0.22 | 1 | 3/4 |
| tnf2_10 | `helix0_E67H` | 91 | 1 | **0.47** | 1/1 | 7.11 | 644 | 0.73 | 0.166 | 1 | 3/4 |
| tnf2_11 | `helix0_E70H` | 91 | 1 | **0.25** | 0/1 | 6.44 | 743 | 0.76 | 0.176 | 1 | 3/4 |
| tnf2_12 | `helix0_R21H` | 91 | 1 | **0.01** | 0/1 | 7.05 | 509 | 0.84 | 0.262 | 1 | 3/4 |

`tnf2_06` / `tnf2_07` are the **matched pair**: identical backbone, switch present and absent.
`tnf2_06` is ranked above its ΔQ because it is the only design with a credible measured
interface — a design that does not bind cannot show pH-dependent release at all.
