# Design specs — the inputs that actually ran

One directory per tool, because the format is the tool's choice and not ours: **BindCraft and
BindCraft2 read JSON settings, BoltzGen and Boltz read YAML entity specs.** Converting either
would break it, so the split is visible in the layout instead of being explained away.

Sixteen files looks like clutter and is not: each is a distinct configuration that ran, and the
differences between them are results. The tables below are the map. Nothing is merged, because
every candidate merge would delete a measurement — the five BindCraft settings share no common
base at all (every field differs between them), and the four early BoltzGen specs are a
deliberate 2×2 that exists so its two factors stay separable.

**Paths are machine-specific.** Every spec carries absolute `/home/ubuntu/...` paths to targets
and output folders on the A100 box the campaigns ran on, kept verbatim: the point of this
directory is to record *what was run*, not to be re-runnable elsewhere without edits. The
target PDBs are not in this repo.

## `bindcraft1/` — the first campaign, and why it was abandoned

Five runs over three axes: the epitope, how much of the target is given, and the budget. This
is the campaign that produced the two failures the whole verifier stack was built to catch —
the top-scoring construct was verbatim ubiquitin, and the top of a 40-design shortlist was
59 % alanine. The submission came from BoltzGen instead.

| | target | hotspots | lengths | final designs | why |
| --- | --- | --- | --- | --- | --- |
| `tnf_smoke.json` | A+B | pocket (6) | 70–100 | 2 | pipeline smoke test, nothing else |
| `tnf_pocket_AB.json` | A+B | pocket (6) | 65–120 | 40 | the standard run |
| `tnf_pocket_ABC.json` | A+B+C | pocket (6) | 65–120 | 40 | same epitope, intact trimer |
| `tnf_lobe_AB.json` | A+B | lobe3 (3) | 65–120 | 40 | the narrowed 3-residue epitope |
| `bc_scale.json` | A+B | pocket (6) | 65–110 | 50 | scale-up of the standard run |

## `bindcraft2/` — the dual-ortholog campaign, all three rounds

Kept in full because the progression between them is the result the experiment card cites:
three rounds, zero accepted designs. JSON has no comments, hence this table.

| | what changed, and why |
| --- | --- |
| `tnf_orthologs_v1.json` | First attempt. Default acceptance filters. Scored **zero** candidates: trajectories died at the screen stage on binder pLDDT 0.59 against a default floor of 0.60. |
| `tnf_orthologs_v2.json` | Relaxed `i_pTM` 0.70 → 0.55 and `i_pAE` 0.35 → 0.45 — affinity proxies, the organisers' lowest-ranked objective — and shortened binders to 65–90 for throughput. Hotspot (≥0.5) and coldspot (≤0.05) gates kept, since those encode epitope focus and cross-reactivity. Still zero: the relaxed *final* filters never saw a candidate, because the deaths were at *screen*. |
| `tnf_orthologs_v3.json` | `target_flexibility` 0.5 + `initial_guess` — rungs 1–2 of BindCraft2's own desperation ladder, which documents these as the measured fix for exactly this failure signature. The ladder cannot self-engage here: it waits 750 trajectories without an acceptance, and one costs ~33 min. `min_plddt_screen` 0.60 → 0.55, since round 1 failed by 0.01; the final gate stays at 0.7 so ProteinMPNN still has to fix the fold. **Designs accepted under these settings are judged against an easier task and are weaker candidates — BindCraft2's own guidance.** |

## `boltzgen/` — the campaign that produced the submission

**The four `screen_*.yaml` are generated**, by [`gen.py screen`](../scripts/gen.py) from its
`CONFIGS` table. They are committed anyway: they are the record of what ran, and the generator
cannot be re-run here because it needs a target PDB this repo does not carry. They screened
three levers — epitope, a helix bias, and length — to lift `design_to_target_iptm` off a
0.11–0.30 plateau. Two chassis survived, `screen_pocket_helix_0` and `screen_pocket_free_1`,
and every design in the submission is built on one of them.

| | epitope | length | secondary structure |
| --- | --- | --- | --- |
| `screen_pocket_free.yaml` | pocket (113/115/143/144) | 60–95 | unconstrained |
| `screen_pocket_helix.yaml` | pocket | 60–95 | 3-helix bundle |
| `screen_pocket5_helix.yaml` | pocket + 67 | 70–110 | 3-helix bundle |
| `screen_lobe3_free.yaml` | lobe3 (88/90/92) | 60–95 | unconstrained |

The remaining four are a **2×2 over method and target**, written that way on purpose: the A+B
dimer truncation depresses interface confidence by ~0.10, so comparing de novo-on-dimer against
nanobody-on-trimer would confound the target with the method. Each cell ran.

| | A+B dimer | intact A+B+C trimer |
| --- | --- | --- |
| **de novo** | `tnf_v2_msa.yaml` | `tnf_denovo_trimer.yaml` |
| **nanobody grafting** | `tnf_nano.yaml` | `tnf_nano_trimer.yaml` |

Two findings live in those files' own comments rather than here, being properties of the specs
themselves: the nanobody track exists because de novo generation collapsed onto ubiquitin, the
strongest small-protein mode in the prior, which is both a non-binder and disqualified by the
de-novo-only filter; and the dimer truncation was **wrong despite preserving the epitope
geometry** — chain C sits 18.7 Å from the hotspots, but TNF-α is only stable as a trimer, so
the folding step returned `target_ptm` ≈ 0.47 for the lone dimer. The model was unsure of the
*target's* own fold, which caps interface ipTM at ~0.15 however good the binder is.
