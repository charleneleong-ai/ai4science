# Per-design evidence

Ranked on **ΔQ floor** (charge gained at the interface going pH 7.4 → 6.0, each histidine
scored on whichever of its three structures gives the *smaller* swing — a lower bound),
then mouse retention, then human affinity.

`robust` = histidines staying inside 6.0–7.6 on **all three** structures; 9 of 26 across the
set do not. `mouse` = interface confidence against mouse TNF-α as a fraction of the human
value, same protocol both sides.

`tnf2_13` filled one of the slots held for a campaign that produced nothing. It is the single
survivor of eight further three-histidine combinations on `helix_0` put through the same
gates — six failed on structure, one on mouse retention at 56 %. It scored 3/4 on novelty, as
predicted from its 2-substitution distance to `tnf2_08`; **all 13 designs pass the ≥3/4 gate.**

| # | design | aa | His | ΔQ floor | robust | mouse | buried Å² | human ipTM | Cys | novelty |
|---|---|---|---|---|---|---|---|---|---|---|
| tnf2_01 | `helix0_R21H_E67H_E70H` | 91 | 3 | **1.38** | 3/3 | 103% | 748 | 0.153 | 1 | 3/4 |
| tnf2_13 | `helix0_R35H_E67H_E70H` | 91 | 3 | **1.35** | 2/3 | 98% | 658 | 0.207 | 1 | 3/4 |
| tnf2_02 | `helix0_R21H_E67H_R72H` | 91 | 3 | **1.14** | 2/3 | 108% | 777 | 0.162 | 1 | 3/4 |
| tnf2_03 | `helix0_E67H_R72H` | 91 | 2 | **1.03** | 1/2 | 88% | 528 | 0.159 | 1 | 3/4 |
| tnf2_04 | `helix0_R21H_E70H` | 91 | 2 | **0.99** | 2/2 | 88% | 889 | 0.175 | 1 | 3/4 |
| tnf2_05 | `helix0_E67H_E70H` | 91 | 2 | **0.89** | 2/2 | 110% | 716 | 0.151 | 1 | 3/4 |
| tnf2_06 | `helix0_R35H_R72H` | 91 | 2 | **0.73** | 1/2 | 74% | 694 | 0.216 | 1 | 3/4 |
| tnf2_07 | `helix0_R72H` | 91 | 1 | **0.57** | 1/1 | 69% | 600 | 0.22 | 1 | 3/4 |
| tnf2_08 | `helix0_E67H` | 91 | 1 | **0.47** | 1/1 | 110% | 644 | 0.166 | 1 | 3/4 |
| tnf2_09 | `free1_T73H` | 74 | 1 | **0.37** | 0/1 | 47% | 923 | 0.634 | 0 | 4/4 |
| tnf2_10 | `free1_parent` | 74 | 0 | **0.0** | n/a | 41% | 802 | 0.497 | 0 | 3/4 |
| tnf2_11 | `helix0_E70H` | 91 | 1 | **0.25** | 0/1 | 82% | 743 | 0.176 | 1 | 3/4 |
| tnf2_12 | `helix0_R21H` | 91 | 1 | **0.01** | 0/1 | 83% | 509 | 0.262 | 1 | 3/4 |

`tnf2_09` / `tnf2_10` are the **matched pair** — identical backbone, switch present and
absent. Both were demoted from the top of the list once the mouse fold showed they retain
only 47 % and 41 % of their human interface, despite the highest conserved-epitope coverage
in the set. They are kept because the pair is still the cleanest control for the switch.
