# Per-design evidence

Every value measured on the design named, not inherited from its parent. pKa is PROPKA
3.5.1 on the mutant in its **bound** context; combinations were measured as combinations.
Cross-reactivity is the fraction of buried interface area *off* the nine human/mouse
divergent positions. Interface confidence is Boltz-2 binder-to-target chain-pair ipTM,
mean across 3 diffusion samples (never best-of-n, never global ipTM).

| # | design | chassis | aa | His | pKa measured in combination | ideal-window His | buried Å² | cross-react | ipTM | clashes | Cys | novelty |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tnf2_01 | `free1_T73H` | free_1 | 74 | 1 | 6.27 | 0 | 923 | 0.79 | 0.634 | 0 | 0 | **4/4** |
| tnf2_02 | `helix0_R21H_E67H_E70H` | helix_0 | 91 | 3 | 7.05 7.05 6.43 | 3 | 748 | 0.64 | 0.153 | 0 | 1 | **3/4** |
| tnf2_03 | `helix0_E67H_R72H` | helix_0 | 91 | 2 | 7.12 6.66 | 2 | 528 | 0.7 | 0.159 | 0 | 1 | **3/4** |
| tnf2_04 | `helix0_R72H` | helix_0 | 91 | 1 | 6.66 | 1 | 600 | 0.78 | 0.22 | 0 | 1 | **3/4** |
| tnf2_05 | `helix0_R35H_R72H` | helix_0 | 91 | 2 | 6.07 6.66 | 1 | 694 | 0.89 | 0.216 | 2 | 1 | **3/4** |
| tnf2_06 | `helix0_R21H` | helix_0 | 91 | 1 | 7.05 | 1 | 509 | 0.84 | 0.262 | 1 | 1 | **3/4** |
| tnf2_07 | `helix0_R21H_E67H_R72H` | helix_0 | 91 | 3 | 7.10 7.12 6.65 | 3 | 777 | 0.8 | 0.162 | 2 | 1 | **3/4** |
| tnf2_08 | `helix0_E70H` | helix_0 | 91 | 1 | 6.44 | 1 | 743 | 0.76 | 0.176 | 1 | 1 | **3/4** |
| tnf2_09 | `helix0_E67H` | helix_0 | 91 | 1 | 7.11 | 1 | 644 | 0.73 | 0.166 | 2 | 1 | **3/4** |
| tnf2_10 | `helix0_E67H_E70H` | helix_0 | 91 | 2 | 7.05 6.44 | 2 | 716 | 0.74 | 0.151 | 2 | 1 | **3/4** |
| tnf2_11 | `helix0_R21H_E70H` | helix_0 | 91 | 2 | 7.04 6.44 | 2 | 889 | 0.7 | 0.175 | 3 | 1 | **3/4** |
| tnf2_12 | `free1_parent` | free_1 | 74 | 0 | — | 0 | 802 | 0.77 | 0.497 | 0 | 0 | **3/4** |

`tnf2_01` and `tnf2_12` are the **matched pair**: identical backbone, switch present and
absent. `tnf2_02`–`tnf2_11` are the titration on the second chassis.

Three candidates were dropped for cross-reactivity below 0.60 despite competitive interface
confidence — including the family's best (0.352 ipTM) which placed 59% of its interface on
divergent residues.
