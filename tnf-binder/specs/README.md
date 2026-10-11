# Design specs — the inputs that actually ran

**Why two formats.** It is not a style choice: BindCraft and BindCraft2 read **JSON** settings
files, Boltz and BoltzGen read **YAML**. Each file is in the format its tool requires, so
converting either would simply break it.

**Paths are machine-specific.** Every spec carries absolute `/home/ubuntu/...` paths to targets
and output folders on the A100 box the campaigns ran on. They are kept verbatim rather than
tidied into relative paths, because the point of this directory is to record *what was run*,
not to be re-runnable elsewhere without edits.

## `bindcraft2/` — the dual-ortholog campaign, all three rounds

The three rounds are kept because the progression between them is the result, and the
experiment card refers to it. JSON has no comments, hence this file.

| | what changed, and why |
| --- | --- |
| `tnf_orthologs_v1.json` | First attempt. Default acceptance filters. Scored **zero** candidates: trajectories died at the screen stage on binder pLDDT 0.59 against a default floor of 0.60. |
| `tnf_orthologs_v2.json` | Relaxed `i_pTM` 0.70 → 0.55 and `i_pAE` 0.35 → 0.45 — affinity proxies, the organisers' lowest-ranked objective — and shortened binders to 65–90 for throughput. Hotspot (≥0.5) and coldspot (≤0.05) gates kept, since those encode epitope focus and cross-reactivity. Still zero: the relaxed *final* filters never saw a candidate, because the deaths were at *screen*. |
| `tnf_orthologs_v3.json` | `target_flexibility` 0.5 + `initial_guess` — rungs 1–2 of BindCraft2's own desperation ladder, which documents these as the measured fix for exactly this failure signature. The ladder cannot self-engage here: it waits 750 trajectories without an acceptance, and one costs ~33 min. `min_plddt_screen` 0.60 → 0.55, since round 1 failed by 0.01; the final gate stays at 0.7 so ProteinMPNN still has to fix the fold. **Designs accepted under these settings are judged against an easier task and are weaker candidates — BindCraft2's own guidance.** |

## Top level — BoltzGen / Boltz

`screen_*.yaml` are the BoltzGen screens that produced the two chassis in the submission;
`tnf_*.yaml` and `tnf_*.json` are target constructions and earlier single-ortholog attempts.
