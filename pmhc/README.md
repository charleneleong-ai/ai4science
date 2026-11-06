# KRAS G12D pMHC binder — Anthropic × Adaptyv Challenge 5

A de novo protein that binds the KRAS G12D peptide `VVVGADGVGK` on HLA-A\*11:01 but not the
wild-type `VVVGAGGVGK`, and not the peptide-free HLA. Specificity outranks affinity. Submissions
close **Sunday 8 November 2026, 23:59 AoE — Monday 9 November, 22:59 Sydney**.

**Status (2026-11-06): all three gates run, the cohort measured twice at five samples, and the
bar calibrated over ten independent runs at +0.106. No design beats the 7OW6 receptor's
peptide-dependence; 4 of 78 are statistically indistinguishable from it and 14 clear zero. The
leader, `mut_proximal_16`, returned +0.112 in both independent cohort runs.** Ranking tier 1 is
not merely unresolvable on this predictor but inverted on the one molecule that can calibrate
it. The replicated, paired wild-type pass then **retracted** the one per-variant result that
had looked mechanistic: it was an artifact of an unpaired difference. The 20-slot submission is
assembled and validated but **not yet uploaded** -- the Proteinbase portal takes entries by
hand -- and nothing is wet-lab measured.

## The finding that set the strategy

The challenge asks for discrimination between two peptides that differ by one residue: Gly12 to
Asp at peptide position 6. That aspartate **points down into the groove**, and its only
solvent-accessible atom is OD1 at 6.4 Å². Casting 262,144 rays from OD1, the fraction that escape
8 Å unobstructed is **0.013%**, or 34 rays — surface area with no approach vector. The probe is
calibrated on the same peptide: position 7's carbonyl, which the T-cell receptor contacts 70
times, returns **39.1%**, and shielded position 4 returns **0.00%**.

The cone depends sharply on how far a side chain must reach, and the dependence is where the
mechanism shows. Over 262,144 rays at a 1.55 Å probe:

| reach | 3 Å | 4 Å | 5 Å | 6 Å | 7 Å | 8 Å |
| --- | --- | --- | --- | --- | --- | --- |
| G12D OD1 | 41.0% | 21.9% | 7.1% | 2.9% | **0.018%** | 0.013% |
| G12D position 6 backbone O | 16.6% | 3.1% | 2.8% | 2.6% | 2.3% | **1.8%** |
| wild-type position 6 backbone O | 18.0% | 5.6% | 3.2% | 3.2% | 3.0% | **2.9%** |

The carboxylate's cone collapses 161-fold between 6 and 7 Å while the backbone carbonyl's holds:
OD1 sits at the bottom of a pocket about 6 Å deep, and the backbone keeps an open channel the
whole way out. An earlier version of this README reported **0.02%** from 4096 rays, which was a
single surviving ray at the sampling floor — the conclusion was right and the number was not.

The one protein known to make this discrimination does not use the side chain either. In 7OW6 the
receptor makes 17 contacts to position 6 within 5.5 Å and **every one is to a backbone atom** —
none to CB, CG, OD1 or OD2.

So the discrimination is indirect readout. Superposing the two complexes on the HLA heavy chain
(0.337 Å residual over shared CAs) leaves the peptide differing only at positions 6–8, by
**0.85 / 0.73 / 0.48 Å** at CA, with the flanks at ≤0.21 Å. Local openness roughly halves:
position 5's carbonyl has a 6.7% escape cone in the mutant against 13.1% in wild type.

Both signals sit below the coordinate error of Boltz-2 and Chai-1. **Selectivity can therefore be
designed toward but not verified**, and is reported as sampled rather than engineered — the same
discipline that kept [`glp1r`](../glp1r) from claiming agonism.

## The other finding: the generator plagiarises the target

The first pilot returned two "designs" that were **β2-microglobulin**, a chain of the target
complex, at 90.5% and 85.7% identity. The target had been truncated to the HLA groove, leaving
β2m's interface vacant, and BoltzGen filled it with β2m.

Those poses scored `design_to_target_iptm` **0.73** and **0.72**, against **0.19–0.29** for the
pilot's genuine designs. Generation moved to the intact complex, which cut the rate from 2 of 2
to **17 of 256**. The challenge requires de novo, zero-shot designs and forbids modifying an
existing binder, so this is a disqualification risk, not a scoring one, and
[`novelty.py`](scripts/novelty.py) is a hard gate run before selection.

**The score bias this was first reported with does not survive the full cohort.** On the pilot's
four designs it looked as though ipTM rewards plagiarism threefold. Over 256 designs it does not:

| group | n | median ipTM | max |
| --- | --- | --- | --- |
| plagiarised (β2m) | 17 | 0.150 | 0.372 |
| copy of another known protein | 16 | 0.186 | 0.880 |
| admissible de novo | 188 | 0.180 | **0.916** |

That is a ratio of **0.83×** — the copies score slightly *worse* — and only **1 of the top 20**
by ipTM is inadmissible. The 3× effect was real but local to target preparation: handed a vacant
natural interface, β2m slots into it and makes a genuinely confident complex; handed the intact
complex there is nothing to slot into. A co-fold score answers "is this a plausible complex", so
it rewards a copy only where the copy makes a better complex than a design would. The gate still
runs before selection — one inadmissible design in the top 20 is one too many when the penalty
is disqualification — but it is a correctness rule, not a bias correction.

## Pipeline

```
256 generated    4 hotspot variants x 64, BoltzGen 0.1.4 against the intact pMHC
239 novelty      < 30% of the design is a verbatim copy of any target or TCR chain
                 17 rejected, every one beta-2 microglobulin, at 30.3-92.6% copied
188 fold         FoldSeek over the designed chain alone, vs CATH50 and the PDB
                 16 more are copies of known proteins the sequence gate could not see
                 35 have an immunoglobulin beta-sandwich as their closest fold
 87 engagement   peptide share >= 0.200, the measured share of the 7OW6 receptor itself
                 78 of those also touch every position their variant nominated
 78 dependence   co-fold screened twice at 5 samples: 318 of 318 predictions scored, 0 failures
                 14 clear zero on the cohort's measured between-run spread (0.018)
                  4 tied with the bar (+0.106, 10 runs), 0 above it
 20 submitted    cap; nothing uploaded yet
```

The denominator stays 256. 188 admissible for 20 slots, so every gate can be conservative for
free.

**The engagement bar is measured, not chosen.** The peptide is only 1.2% of the complex's
accessible surface, so gripping the HLA helices and ignoring the peptide is the easy thing for a
generator to do — and it fails ranking tier 2, peptide-dependence. The 7OW6 receptor says what
engagement actually looks like here: **55 contacts on the peptide against 220 on the HLA, a share
of 0.200**, with just 5 contacts to the mutated position 6. A plausible-sounding gate of "half
the interface on the peptide" would have rejected the one binder known to work.

Scored against that bar the four variants separate sharply, which is the strongest argument for
having run all four:

| variant | nominated positions | median peptide share | admissible | engaged |
| --- | --- | --- | --- | --- |
| `mut_proximal` | 5, 6 | **0.206** | 48 | 28 |
| `bulge_wide` | 4, 5, 6, 7, 8 | 0.161 | 48 | 23 |
| `bulge_core` | 5, 6, 7 | 0.046 | 50 | 19 |
| `hub` | 7, 8 | 0.029 | 42 | 17 |

Nominating the mutation and its neighbour is what produces peptide-reading designs; nominating
the hub at 7–8 produces groove-grippers, even though position 7 is the most exposed residue in
the epitope and the receptor's own contact hub.

An earlier version of this table reported 0.207 / 0.167 / 0.051 / 0.032, taking the upper of the
two middle values of each 64-design group rather than their mean. The ordering and the
conclusion are unchanged; `mut_proximal` still clears the 0.200 bar at its median.

## The control says tier 1 is unreachable, and now says it with a number

The co-fold screen put the 7OW6 receptor through all three conditions. It is the only molecule
known to discriminate these two peptides, and it does not discriminate them here:

| 7OW6 receptor, Boltz-2 | ipTM |
| --- | --- |
| G12D complex | 0.957 |
| **wild-type complex** | **0.958** |
| peptide-free HLA | 0.842 |

That difference is not merely too small to use -- it has the wrong sign, and repetition makes
the sign solid rather than softening it. Paired within each of eight independent runs, the
control scores wild type **above** G12D by +0.0027 +- 0.0006, in seven of the eight. The one
receptor known to discriminate these peptides is scored here as preferring the peptide it does
not recognise. The sub-angstrom backbone shift predicted that tier 1 would be unreachable; the
control shows the scoring function is not just blind to it but inverted on it, so no
selectivity a design appears to show on this predictor can be read as signal.

Peptide-dependence does resolve, and its bar is now measured rather than read off one run.
Ten independent invocations of the control put the drop at **+0.106**. Eight of them, run as
small batches, span only +0.095 to +0.108 with a spread of 0.004 against the 0.007 each run
reports internally -- a permutation test over the pooled samples places that spread at the 4th
percentile of the no-run-effect null, so for the control a single run's own standard error is
already honest. The spread of 0.016 over all ten is one run's doing, diagnosed below. Against
that bar, the two five-sample cohort runs:

| | peptide-dependence drop |
| --- | --- |
| `mut_proximal_16` | **+0.112** (+0.112, +0.112) |
| `mut_proximal_26` | **+0.095** (+0.109, +0.081) |
| the 7OW6 receptor, 10 runs | +0.106, spread 0.016 |
| best of the rest (`hub_62`) | -0.016 +- 0.023 |
| median of 78 | **-0.010** |
| worst | -0.121 |

**14 of 78** clear their own standard error at 2 sigma. The median design scores slightly
*better* against the peptide-free receptor than against the complex, which is what a
groove-gripper looks like when the peptide it ignores is taken away. Both leaders come from
`mut_proximal`, the variant the engagement gate already favoured on the generator's own poses
-- two independent measurements agreeing about which nomination strategy worked.

Whether either leader beats the control is now settled, and the answer is no. Replicating the
whole 159-spec pass gives every design two independent measurements, and the field's
reproducibility is the denominator each one has to borrow: a between-run SD of **0.018** against
a within-run error of 0.015, leaving a run effect of 0.009. Two runs cannot estimate a single
design's own spread -- the estimate carries ~70% relative error, and a pair that happens to
agree reports zero -- so every verdict divides by the pooled number instead. On it: **0 designs
above the bar, 4 tied, 74 below**, and 14 of 78 clear zero. The verdict is the same whether the
bar is the 10-run +0.106 or the +0.101 that excludes the degraded-pose run.

The replicate also ran the artifact test, and the leaders split on it. `mut_proximal_16`
returned **+0.112 in both runs** -- a between-run SD of 0.0001, so its peptide-dependence is not
a lucky apo pose. `mut_proximal_26` fell from +0.109 to +0.081, which is the bias doing exactly
what was predicted of it.

**The bar moved once, and the reason matters more than the number.** The 159-spec cohort run
scored the control's drop at +0.150, four of its own standard errors above the probe's +0.104.
It is not run-to-run noise: the spec is byte-identical, the processed structure hashes equal,
the flags match, and eight replicates span only +0.095 to +0.108. What degraded was the
*peptide-free* pose alone -- pTM 0.836 against 0.877, and receptor-to-HLA pair ipTM 0.710
against 0.786, while pLDDT barely moved and the complex condition reproduced to three decimals.
With no peptide in the groove the receptor's placement is underdetermined, which is why the apo
term's spread is 4-10x the complex term's, and why the drop inherits all its noise from the
condition that defines it. The 12 designs measured in both runs shifted their apo scores by
+0.0006 +- 0.0150 while the control's fell 0.045, so this was not a run-level offset that an
internal control would absorb -- which is the reason the cohort's +0.150 is not used as the bar.

One consequence is a bias, and it runs toward false positives: an unlucky apo pose *inflates* a
design's drop. The top of this ranking is therefore enriched for peptide-free-pose failures as
well as for genuine peptide reading, and the replicate pass is the test -- artifacts regress
toward the middle, real peptide-readers do not.

**Reproducing the ranking.** `pool-replicates` writes the table of record and `submission.py`
reads it:

```bash
python scripts/cofold_followup.py pool-replicates \
  --out-dir cofold_ctrl_out --out-dir cofold_full5_out --out-dir cofold_full5b_out \
  --out cofold_spread_pooled.csv
python scripts/submission.py build          # defaults to that table
```

**Ranking tier 1.** `select-wildtype` writes both peptides for the finalists *and* for the
control, and the pass is run at least twice as separate invocations:

```bash
python scripts/cofold_followup.py select-wildtype --out cofold_wt_specs
bash scripts/run_cofold.sh cofold_wt_specs cofold_wt_out 5 2
python scripts/cofold_followup.py report-wildtype --out-dir cofold_wt_out \
  --out cofold_wildtype.csv
```

Both conditions in one pass because the difference is only interpretable paired: an invocation
shifts both conditions together, so forming the difference inside it cancels that shift, while
taking the two sides from different invocations does not — and the error bar cannot see what it
left out. The control travels with the finalists for the same reason the bar does everywhere
else here: it is the threshold, and it has to be measured in the same runs.

**What it returned, and what it retracted.** 420 predictions, 0 failures. The control's own
difference is **+0.0018 ± 0.0065** over the two invocations — still positive, so the predictor
still prefers the peptide the receptor does not recognise, but tied zero at this sample size.
3 of 20 finalists sit below that bar, 13 tied, 4 above; the cohort averages +0.0090 ± 0.0058.
An earlier *unpaired* pass over the same designs had shown `bulge_core` — the variant placed
over the bulge — 4 of 4 negative at mean −0.0149, Fisher p = 0.0072, which was the one result
here that looked mechanistic. Paired, those four average **+0.0115** at **p = 0.1011**, because
`bulge_core_27` moves across zero. The measured order is `mut_proximal` < `bulge_wide` <
`bulge_core` < `hub`, not the pre-registered one; only `hub` being the most positive group
survives. Reported, never gated, and no design was re-ranked on it.

**Rebuilding the deliverables.** Both render from the committed files, so neither can carry a
stale number:

```bash
uv run --extra artifact python scripts/build_artifact.py   # -> docs/artifact/index.html
cd docs && pdflatex report.tex && pdflatex report.tex      # -> docs/report.pdf
```

The W&B side is the same files again, logged rather than rendered. The ledger goes up as one run
per experiment, and the report's Results and run-log sections are rewritten in place from it,
leaving the narrative sections alone:

```bash
uv run --extra wandb python scripts/log_wandb.py --entity <entity>
uv run --extra wandb python scripts/build_wandb_report.py create --entity <entity>
uv run --extra wandb python scripts/build_wandb_report.py refresh --entity <entity> --url <report url>
```

Every row names the scale its error bar is built on, in a `sem_basis` column, because the same
`drop_sd` header means a design's own sample spread in a single-pass table and the cohort's
between-run spread in the pooled one. The published evidence file carries it, so the ranking's
denominator is readable without knowing which command wrote the file.

The ranking within that list is not claimed from the screen. Each of its drops is a difference
of two single diffusion draws, and a 26-design probe at five samples settled the question the
wrong way: two of the screen's top four **reverse sign** when measured five times. The whole
cohort is therefore re-run at five samples and read by
[`cofold_followup.py`](scripts/cofold_followup.py) `report-spread`, which owns the campaign's one
significance rule — a difference must clear two standard errors of the mean, and `drop_sem`, not
the per-draw `drop_sd`, is that denominator. Designs inside the band are reported as ties.

The fold gate cannot be a TM-score cutoff. **All 256 designs match something at TM ≥ 0.5 bar
two**, because an 80-residue helical bundle resembles one of the thousands already deposited, so
a TM ≥ 0.5 reject would discard the cohort. What separates them is the sequence identity to that
structural match: **85–91%** to a PDB chain is a known protein handed back, while 0–16% at the
same TM is a de novo bundle sharing a common fold. The gate rejects a high-TM match at ≥30%
identity, plus the target's own immunoglobulin β-sandwich class (CATH 2.60.40), which is β2m's
architecture and the class an earlier submission of ours failed the organisers' check on at
**37.5%** identity while a helical design passed at **38.1%**.

That 30% line is why the gate matters rather than merely agreeing with the sequence gate. The 16
designs it catches alone sit at **TM 0.93–0.98 and 30–43% identity** to their structural match —
the exact band those two recorded outcomes straddle — and several match CATH 2.60.40.10 directly.
The sequence gate cannot see them because it compares only against the six target-derived
chains, and these designs resemble β2m as deposited in *other* entries.

## Layout

| Path | What it does |
| --- | --- |
| [`scripts/target_sites.py`](scripts/target_sites.py) | carves the TCR off both depositions and measures per-residue peptide exposure |
| [`scripts/epitope_geometry.py`](scripts/epitope_geometry.py) | surface split, per-atom SASA, and how the known receptor reads position 6 |
| [`scripts/compare_conformers.py`](scripts/compare_conformers.py) | superposes the two complexes and localises the mutation's backbone effect |
| [`scripts/reachability.py`](scripts/reachability.py) | the escape-cone probe behind the 0.013% result |
| [`scripts/prepare_target.py`](scripts/prepare_target.py) | writes the generation target, the three co-fold conditions, novelty references and specs |
| [`scripts/novelty.py`](scripts/novelty.py) | the hard first gate against regurgitated target chains |
| [`scripts/engagement.py`](scripts/engagement.py) | peptide-vs-groove contact split, with the 7OW6 receptor as the bar |
| [`scripts/cofold_specs.py`](scripts/cofold_specs.py) | applies all three gates, records every design's fate, writes the co-fold specs |
| [`scripts/cofold_scores.py`](scripts/cofold_scores.py) | reads the screen and reports each design's peptide-dependence drop |
| [`scripts/run_fold_gate.sh`](scripts/run_fold_gate.sh) | the detached fold-gate run over the whole cohort |
| [`scripts/run_cofold.sh`](scripts/run_cofold.sh) | the detached co-fold screen; a fourth argument repeats the pass as separate invocations |
| [`data/cofold_spread_pooled.csv`](data/cofold_spread_pooled.csv) | **the table the submission ranks on** — both cohort runs pooled, `drop_sem` on the cohort's between-run spread |
| [`data/cofold_replicates.csv`](data/cofold_replicates.csv) | one row per design per invocation, the audit trail under the pooled table |
| [`docs/experiments/pmhc/verdicts.yaml`](docs/experiments/pmhc/verdicts.yaml) | one row per experiment: gate, result, verdict, and what is not established |
| [`scripts/cofold_followup.py`](scripts/cofold_followup.py) | the five-sample cohort pass, the one significance rule the campaign uses, the replicated bar and inverted tier-1 measurement, and the wild-type finalists |
| [`scripts/submission.py`](scripts/submission.py) | ranks on the five-sample means and writes the upload CSV, refusing to rank without the control |
| [`scripts/fold_novelty.py`](scripts/fold_novelty.py) | FoldSeek fold gate, for the copies sequence identity cannot see |
| [`tests/test_novelty.py`](tests/test_novelty.py) | the admissibility gates, against copies they have to catch |
| [`scripts/run_generate.sh`](scripts/run_generate.sh) | the detached BoltzGen sweep |
| [`docs/experiment-card.md`](docs/experiment-card.md) | frozen inputs, versions, budget, and what is not established |
| [`docs/experiments/pmhc/campaign.yaml`](docs/experiments/pmhc/campaign.yaml) | the measured numbers that are not a CSV row: the escape-cone grid, the funnel, the bar, the variance components |
| [`scripts/build_artifact.py`](scripts/build_artifact.py) | renders the progress page from those files; nothing on the page is typed in |
| [`scripts/log_wandb.py`](scripts/log_wandb.py) | mirrors the ledger to W&B, one run per experiment, so the runs table *is* the verdict table |
| [`scripts/build_wandb_report.py`](scripts/build_wandb_report.py) | the project's one W&B report: narrative written once, Results and run log regenerated in place |
| [`docs/artifact/index.html`](docs/artifact/index.html) | the built page, committed on every rebuild |
| [`docs/report.tex`](docs/report.tex) · [`docs/report.pdf`](docs/report.pdf) | the written report: method, funnel with every denominator, and what is not established |

## Running it

Measurements, on a host with the BoltzGen venv (`gemmi`, `biopython`):

```bash
cd ~/pmhc
for id in 7OW6 7OW5 7OW4; do curl -sSLo data/$id.cif https://files.rcsb.org/download/$id.cif; done
export PYTHONPATH=$HOME/pmhc
python scripts/target_sites.py --data data        # writes the carved *_pmhc.cif files first
python scripts/epitope_geometry.py --data data
python scripts/compare_conformers.py --data data
python scripts/reachability.py --data data
python scripts/prepare_target.py --data data
```

`target_sites.py` must run first: it writes the carved `7OW6_pmhc.cif` and `7OW5_pmhc.cif` that
the other three read.

Generation, detached, with the launch gate and a heartbeat:

```bash
setsid nohup bash scripts/run_generate.sh </dev/null >>logs/gen_$(date -u +%Y%m%dT%H%M%SZ).log 2>&1 & disown
setsid nohup autoresearch-heartbeat --label pmhc-gen --pattern '[r]un_generate.sh' \
  --log "$(ls -t logs/gen_*.log | head -1)" --interval 1800 --stall-after 3600 \
  --state-cmd 'find $HOME/pmhc/gen_out -name "*.cif" | wc -l' </dev/null >>logs/heartbeat.log 2>&1 & disown
```

Then the novelty gate on each variant's metrics:

```bash
python scripts/novelty.py gen_out/bulge_core/final_ranked_designs/all_designs_metrics.csv
```

## What is not established

- No wet-lab measurement: no K<sub>D</sub>, no expression, no sensorgram.
- Selectivity against the wild-type peptide is not verifiable with these tools, for the reason
  above. It is sampled, not engineered.
- The aspartate rotamer every design is built against comes from a receptor-bound structure.
  7OW4, the receptor-free G12D crystal, has no coordinates for positions 5 and 6 at all.
- The screen runs one diffusion sample per condition, so the reported drop is a difference of
  two single draws. Until [`cofold_followup.py`](scripts/cofold_followup.py) has sized that
  spread, no ranking by drop is claimed beyond clearing the control's.
- `design_peptide_pair_iptm` is reported for audit and is not a gate: on the control the pair
  value reaches 0.894 between two chains that make no contact at all.
- No ranking between surviving designs is claimed until both verifiers have run and their
  agreement is measured.
