# Experiment card — KRAS G12D pMHC binder (Challenge 5)

Frozen 2026-11-02. Generation in flight; every number below is measured, and rows that are not
yet measured say so rather than carrying a placeholder.

## Task

Design a protein that binds the KRAS G12D peptide `VVVGADGVGK` presented on HLA-A\*11:01, and
does **not** bind either the wild-type peptide `VVVGAGGVGK` on the same allele or the same HLA
with an empty groove. Track 3, 20 designs, 10–250 aa, de novo and zero-shot only. Closes
Sunday 8 November 23:59 AoE, which is Monday 9 November 22:59 in Sydney.

Ranking is tiered, specificity first: mutation specificity, then peptide-dependent binding, then
affinity (lowest K<sub>D</sub>). The readout is a sensorgram per condition; binding is a fittable
K<sub>D</sub>, or an association shift above 300% of the negative control.

## Inputs, frozen

| Input | Value |
| --- | --- |
| Mutant complex | PDB [7OW6](https://www.rcsb.org/structure/7OW6), chains A/B/C, TCR chains D/E deleted |
| Wild-type complex | PDB [7OW5](https://www.rcsb.org/structure/7OW5), same treatment |
| Receptor-free mutant | PDB [7OW4](https://www.rcsb.org/structure/7OW4) — **peptide positions 5 and 6 have no coordinates** |
| HLA heavy chain | 7OW6 chain A, 276 residues as deposited |
| β2-microglobulin | 7OW6 chain B, 99 residues |
| Peptide | chain C, 10 residues, position 6 Asp (mutant) / Gly (wild type) |
| Source paper | Poole et al., *Nat Commun* 2022 |

One heavy chain and one β2m serve all three co-fold conditions, which differ only in the peptide.
The two depositions disagree by a leading construct residue (276 vs 275); carrying that into the
comparison would have placed a sequence difference inside the exact contrast selectivity is read
from, so `prepare_target.py` uses the mutant's chains throughout and asserts the two alleles
differ only at the construct ends.

## Tools, exact versions

| Tool | Version | Role |
| --- | --- | --- |
| BoltzGen | 0.1.4 | generator, `--protocol protein-anything` |
| Boltz | 2.2.1 | co-fold verifier (independent of the generator's own folding) |
| Chai-1 | `chai.sif`, built 2026-10-16 | second co-fold verifier |
| gemmi | 0.6.5 | structure I/O, superposition |
| biopython | 1.88 | Shrake–Rupley SASA, pairwise alignment |
| PROPKA / pKAI | 3.5.1 / 1.2.0 | available, not yet used on this target |
| torch | 2.14.1 | A100 80 GB, single device |

## What was measured before generating

The decisive measurements, all on the carved complexes, with the scripts that produced them:

| Measurement | Value | Script |
| --- | --- | --- |
| Peptide share of complex accessible surface | 239.5 Å² of 19,612 Å² — **1.2%** | `epitope_geometry.py` |
| Position-6 side-chain SASA (mutant) | 8.8 Å²; OD1 6.43, OD2 0.00, CG 1.21 | `target_sites.py` |
| Escape cone over OD1, 8 Å reach, 1.55 Å probe | **0.013%** — 34 of 262,144 rays | `reachability.py` |
| Escape cone over OD1 at 3 / 4 / 5 / 6 / 7 Å | 41.0 / 21.9 / 7.1 / 2.9 / 0.018% | `reachability.py` |
| Escape cone, position-6 backbone O, 8 Å (mutant / wild type) | 1.82% / 2.86% | `reachability.py` |
| Escape cone, position-7 carbonyl (probe control) | 39.1% | `reachability.py` |
| Escape cone, position-4 carbonyl (probe control) | 0.00% | `reachability.py` |
| TCR contacts to position 6 at 5.5 Å | 17, **all to backbone**, none to CB/CG/OD1/OD2 | `epitope_geometry.py` |
| TCR interface split at 4.5 Å (the engagement bar) | 55 peptide / 220 HLA — a peptide share of **0.200** | `engagement.py --control` |
| TCR contacts to position 7 at 5.5 Å | 70 of 182 total | `epitope_geometry.py` |
| Peptide backbone shift, mutant vs wild type | CA 0.85/0.73/0.48 Å at positions 6/7/8; flanks ≤0.21 Å | `compare_conformers.py` |
| HLA superposition residual (the baseline that shift sits against) | 0.337 Å | `compare_conformers.py` |

These were taken **before** any design was generated, and they set the strategy rather than
rationalising it afterwards.

One correction to how they were first reported. The escape cone was initially quoted as
**0.02% of 4096 rays**, which is a single surviving ray — the sampling floor, not a measurement,
and the most pessimistic corner of a probe × reach grid that had not been run. The grid has since
been run at 262,144 rays: the value is 0.013% (34 rays), the result is robust across probes
1.40–2.00 Å, and the reach dependence has a 161-fold cliff between 6 and 7 Å that locates the
carboxylate at the bottom of a pocket about 6 Å deep. The conclusion did not change; the evidence
for it did. `reachability.py` now reports `free_rays` beside the fraction, so a number resting on
one ray cannot be quoted again without that being visible in the table.

## Sampling budget

Generation is 4 hotspot variants × 64 designs = **256 attempted**, against the intact complex.
The denominator stays 256 in every funnel reported, including the designs rejected for not being
designs at all.

| Variant | Hotspots (peptide positions) |
| --- | --- |
| `bulge_core` | 5, 6, 7 |
| `bulge_wide` | 4, 5, 6, 7, 8 |
| `mut_proximal` | 5, 6 |
| `hub` | 7, 8 |

An earlier pilot of 2 designs against a **truncated** target (α1+α2 and the peptide, β2m removed)
is excluded from that budget and reported separately below, because the target differed.

## Filters, and when they were fixed

| Gate | Criterion | Fixed before or after seeing results |
| --- | --- | --- |
| novelty | < 30% of the design is a verbatim copy of any target or TCR chain | **After** the pilot — see the incident below. Declared, not pre-registered. |
| fold novelty | top FoldSeek hit must not be a known protein at high identity, nor the target's own immunoglobulin β-sandwich class | **After** the smoke batch. Declared, not pre-registered. |
| engagement | peptide share ≥ 0.200 and every nominated hotspot touched, on the generator's own pose | Criterion **before**; the 0.200 figure is the 7OW6 receptor measured with the same code, so the threshold is calibrated rather than chosen |
| peptide-dependence | co-fold interface against `apo_hla` must degrade relative to `complex_g12d` | Before |
| selectivity | co-fold against `complex_wt` | **Reported, not gated.** Not a judgement call: on the one molecule that can calibrate it the predictor's sign is inverted — see below |
| ranking denominator | the cohort's pooled between-run SD, not a design's own | **After** the cohort was replicated. Declared, not pre-registered: two runs per design cannot estimate a per-design spread |

## The novelty incident

The pilot's first two designs were β2-microglobulin, a chain of the target, at **90.5% and 85.7%
identity**. The generation target had been truncated to the HLA groove, leaving β2m's own
interface vacant, and the generator filled it with the protein that belongs there. Those poses
scored `design_to_target_iptm` **0.73 and 0.72**; genuinely novel designs on this target score
**0.16–0.29**.

Two consequences, both load-bearing:

1. Generation moved to the intact complex. This reduced the rate from 2 of 2 to **17 of 256**,
   and did not eliminate it, so novelty is a hard gate rather than a precaution.
2. **The novelty gate runs before selection**, because admissibility is binary and the penalty is
   disqualification rather than a worse score.

### The score bias was overstated, and the cohort withdraws it

Point 2 was first justified by a stronger claim — that ipTM rewards plagiarism threefold, so a
score-first ranking selects inadmissible designs preferentially. That was measured on the pilot's
four designs. Over all 256 it does not hold:

| group | n | median ipTM | mean | max |
| --- | --- | --- | --- | --- |
| plagiarised, sequence gate rejects (all β2m) | 17 | 0.150 | 0.168 | 0.372 |
| copy of another known protein, fold gate only | 16 | 0.186 | 0.285 | 0.880 |
| immunoglobulin β-sandwich fold | 35 | 0.163 | 0.204 | 0.481 |
| admissible de novo designs | 188 | 0.180 | 0.280 | **0.916** |

The ratio of medians is **0.83×**: the β2m copies score slightly *worse* than admissible designs,
39 of 188 admissible designs beat the best plagiarised one, and **1 of the top 20** by ipTM is
inadmissible (`mut_proximal_31`, caught by the fold gate at TM 0.980 / 35.1% identity to
`7v6q`). The earlier "genuine designs score 0.16–0.29" described two pilot designs, not the
design space, which reaches 0.916.

The 3× effect was real but **local to target preparation**: given a vacant natural interface the
copy slots into it and makes a genuinely confident complex, and given the intact complex there is
nothing to slot into. A co-fold score answers "is this a plausible complex", never "did you
design this", so it favours a copy only where the copy would make a better complex than a design.
Nothing about the pipeline changes — novelty still gates first — but the justification is
correctness, not bias.

The 30% threshold was chosen after seeing 86–91% and 0% as the observed two clusters. It is
declared here as post-hoc, and no statistic derived from it is presented as pre-registered.

### The gate did not measure what it claimed

The gate shipped as a **global** alignment's matches over the shorter sequence, and that cannot
represent a truncated copy at all — the exact case its own docstring named. Measured on verbatim
slices of 7OW6: a 40-residue slice of the 276-residue heavy chain scored **0.000**, and a
60-residue slice of the TCR β chain, the known receptor the rules forbid modifying, also scored
**0.000**. Both would have been written out as novel. How badly it underestimates scales with the
ratio of the lengths, which is why the hole is invisible on a shortened test reference.

The gate is now the fraction of the design that is a verbatim copy of some target chain, by local
alignment, normalised by the design's own length. On controls the same slices score 1.000. The
replacement was chosen on **synthetic controls alone**, independent of the cohort, and the
threshold stayed at the 30% already declared — but the metric change is post-hoc and is declared
as such here. The superseded number is kept in the output as `whole_chain_identity` so the change
is auditable, and [`tests/test_novelty.py`](../tests/test_novelty.py) pins both the new behaviour
and the old blindness against the real chains at their deposited lengths.

On the 192 designs generated so far the two metrics do not agree: the copied fraction rejects
**10**, the superseded one rejects 15, and its five extra rejections score 0.074–0.256 copied.
Its genuine-design p95 is 0.435 against a 0.30 threshold, so it was operating inside its own
noise band; the copied fraction separates cleanly, with genuine designs at ≤0.256 and every
rejection at ≥0.303.

### The fold gate, and a rule that was measured and withdrawn

Sequence identity is necessary and not sufficient: both smoke designs that cleared it at 0%
identity were same-fold matches at **TM 0.74–0.78**. But a TM cutoff cannot be the gate either —
all 11 designs of the first batch matched something at TM 0.71–0.93.

Run over the whole cohort the gate is not redundant with the sequence gate: it independently
rejects all 17 designs the sequence gate caught, **and 16 more** at TM 0.93–0.98 with 30–43%
identity to their structural match. The sequence gate cannot see those, because it compares only
against the six target-derived chains while these designs resemble β2m as deposited in other
entries. A further **35** survivors have an immunoglobulin β-sandwich as their closest fold.
Funnel: **256 generated → 239 after sequence → 188 after fold.**

The 30% identity line for `copy_of_known` was set by analogy to the sequence gate's threshold and
is justified by evidence **external** to this cohort — the organisers' own two recorded outcomes
at 37.5% (failed) and 38.1% (passed). It is not tuned on these designs, but it is not
pre-registered either, and it rejects 16 designs that a more permissive line would keep. With 188
survivors for 20 slots that caution is free, which is the only reason it is set where it is.

One rule was added and then withdrawn on measurement. Ranking each design's hits by how
disqualifying they are, rather than by TM, so that a lower-TM immunoglobulin match would outrank
a strong helical one, labelled **8 of 11** designs immunoglobulin β-sandwich at TM 0.42–0.72 —
including one whose own best matches were TM 0.865 and 0.995 to other folds. The search runs at
`-e 10000` over 1,000 hits per query by design, so some 2.60.40 entry appears for nearly every
query; reading the tail of that list reports the most spurious hit, not the most informative. The
top-TM hit defines the fold, and the class check applies to it.

## The measurement of record

The one-sample screen is retrospectively a pilot. Its drops were differences of two single
diffusion draws, and a 26-design probe at five samples **reversed the sign** of two of its top
four. Everything below is five diffusion samples per condition, and the whole 159-spec cohort
pass was run **twice independently** — 318 predictions, 318 scored, 0 failures.

### The bar, calibrated over ten independent invocations

| | drop |
| --- | --- |
| control, 8 small-batch invocations | +0.101, spread 0.004 (span +0.095 to +0.108) |
| control, 26-spec probe | +0.104 |
| control, 159-spec cohort run 1 | **+0.150** |
| control, 159-spec cohort run 2 | +0.101 |
| **pooled over all 10** | **+0.106**, spread 0.016 |

Nine of ten measurements of the same molecule agree; one does not. That run is excluded from
the +0.101 figure on a **structural** criterion measured independently of its drop, not on the
drop itself: its peptide-free prediction scored pTM 0.836 against 0.877 and receptor-to-HLA pair
ipTM 0.710 against 0.786, while pLDDT barely moved and the complex condition reproduced to three
decimals. It is a docking failure in one condition of one run.

An internal control does not excuse it, which was measured rather than assumed: across the 12
designs present in both the probe and cohort runs, design apo scores moved **+0.0006 ± 0.0150**
while the control's fell **0.045**. The shift was not shared, so pairing designs against that
run's control would have charged every design for one bad prediction.

### Two variance components, and the one that is not there

`Var(run means)` already contains the within-run sampling error, so the two **subtract**:
`run_effect = sqrt(between_sd² − within_sem²)`. On the control, between 0.004 against within
0.007 gives a negative variance, and a 20,000-draw permutation test over the pooled samples puts
the observed spread at the **4th percentile** of the no-run-effect null. There is no resolved
between-run component for the control, and its single-run standard error is if anything
conservative. An earlier version of this analysis combined the two in quadrature; that
double-counted and is withdrawn.

Designs are different: pooled between-run SD **0.018** against a within-run error of 0.015,
leaving a run effect of **0.009**. Every per-design verdict divides by the pooled figure over
√runs, because two runs cannot estimate one design's own spread — the estimate carries ~70%
relative error and a pair that happens to agree reports 0.000, which would make a +0.011 drop
read as significant and a reproducible +0.112 read as infinitely so. The control is excluded
from the pool it is the bar for, so the bar cannot widen its own margin.

### Ranking tier 1 is inverted, not merely unresolved

Paired within each invocation, the control scores **wild type above G12D by +0.0027 ± 0.0006**
over the eight invocations that predicted `complex_wt`, positive in seven of the eight. (The
bar's ten invocations include the two cohort runs, which carried no wild-type condition, so tier
1 rests on eight and not on ten.) The receptor known to discriminate these peptides is scored as
preferring the peptide it does not recognise. The
sub-ångström backbone shift predicted that tier 1 would be unreachable; this says the scoring
function is not blind to it but anti-correlated with it, so no selectivity a design shows in
this column can be read as signal in either direction.

### The result

| | peptide-dependence drop |
| --- | --- |
| `mut_proximal_16` | **+0.112** (+0.112, +0.112 — between-run SD 0.0001) |
| `mut_proximal_26` | +0.095 (+0.109, +0.081) |
| the 7OW6 receptor, 10 runs | +0.106 |
| median of 78 | −0.010 |
| worst | −0.093 |

**0 designs above the bar, 4 tied, 74 below. 14 of 78 clear zero.** Unchanged whether the bar
is the 10-run +0.106 or the +0.101 that excludes the degraded-pose run, which is the most that
can be said for it. All 20 slots are used: the cap is 20, the admissible field is 78, and a
submission of 20 with 4 tied and 0 above is reported as exactly that.

### The finalists against the wild-type peptide

All 20 submitted designs were co-folded against `VVVGAGGVGK`. **Both peptide conditions are
written into a single invocation** for every finalist and for the control, and the whole pass was
run twice as independent invocations — 420 predictions, 0 failures. The pairing is the point: the
control's bar is formed inside one invocation, so a design's difference has to be formed the same
way or the run offset does not cancel out of it.

The control's own difference is **+0.0018 ± 0.0065 over the 2 invocations** (+0.0030, +0.0007).
Still positive — the predictor prefers the peptide the receptor does not recognise — but tied
zero at this sample size, so the bar is weak and so is any ranking against it. Of 20 finalists,
3 sit below it, 13 are tied and 4 are above; 7 of 20 are negative in raw sign. The cohort
averages **+0.0090 ± 0.0058**, a spread over designs rather than over runs. Between-run SD is
0.0092 against a within-run 0.0093, leaving a run effect of ~0.

**The earlier, unpaired pass's headline did not reproduce, and that is the finding.** Its
difference took the wild-type mean from one invocation and the G12D mean from others. Measured
that way, `bulge_core` — the variant whose interface was placed over peptide positions 5–8, the
bulge itself — looked 4 of 4 negative at mean −0.0149, Fisher exact p = 0.0072 against 3 of 16
elsewhere. Paired, the same four designs average **+0.0115** at **p = 0.1011**, because
`bulge_core_27` moves across zero to +0.0767.

| variant | hotspots | n | mean wt − G12D (paired) | unpaired claim |
| --- | --- | --- | --- | --- |
| `mut_proximal` | 5, 6 | 5 | −0.0012 | −0.0073 |
| `bulge_wide` | 4, 5, 6, 7, 8 | 8 | +0.0106 | +0.0052 |
| `bulge_core` | 5, 6, 7 | 4 | **+0.0115** | **−0.0149** |
| `hub` | 7, 8 | 3 | +0.0183 | +0.0209 |

The measured order is `mut_proximal` < `bulge_wide` < `bulge_core` < `hub`, not the
pre-registered one. The only ordinal fact that survives is that `hub`, which puts the least area
on the peptide, is the most positive group (0 of 3 negative).

**This is not evidence of selectivity, and no design was re-ranked on it.** The bar's sign is
inverted on the only molecule known to discriminate these peptides; the bar is tied zero on two
invocations; and the variant grouping was chosen after reading the per-design table, so it is
post-hoc either way. What the replication bought was the retraction of a result that looked
mechanistic and was an artifact of its own denominator.

### The bias the replicate was run to test

An unlucky peptide-free pose *inflates* a drop, because the drop is complex minus apo. The top
of this ranking is therefore enriched for apo-prediction failures as well as for genuine peptide
reading, and the apo term carries 4–10× the complex term's spread precisely because removing the
peptide leaves the binder's placement underdetermined. Replication is the test: artifacts
regress toward the middle, real signal does not. The two leaders split on it — `mut_proximal_16`
returned +0.112 twice, `mut_proximal_26` fell from +0.109 to +0.081.

## What is not established

- **Selectivity is not engineered, and cannot be verified with these tools.** The only atom that
  chemically distinguishes the two targets is OD1, and it has a 0.013% escape cone — surface area
  with no approach vector. The discrimination available is indirect: a 0.5–0.85 Å backbone shift
  and roughly halved local openness at positions 5–6. Both are below the coordinate error of
  Boltz-2 and Chai-1, so a co-fold comparison between `complex_g12d` and `complex_wt` cannot
  resolve them. Selectivity is therefore sampled and reported, never claimed as designed.
- **The aspartate rotamer any design must read has never been observed unliganded.** 7OW4, the
  receptor-free G12D crystal, has no coordinates for positions 5 and 6; 7OW6 holds the peptide in
  a conformation a TCR selected. Every design is built against a receptor-shaped conformation.
- **No wet-lab measurement of anything here.** No K<sub>D</sub>, no expression, no sensorgram.
- **Affinity is not calibrated, and tier 3 is not addressed at all.** The quantity measured here
  is peptide-*dependence*, not affinity. No K<sub>D</sub> is predicted and nothing ranks designs
  by binding strength; a design tied with the control on dependence may bind far more weakly.
- **Only one co-fold predictor was used.** The standing discipline on this project is that every
  load-bearing claim gets an independent second method, and that is **not met here**: the Chai-1
  container was never located, so the consensus pass did not run. Every number above is Boltz
  2.2.1 scoring structures Boltz 2.2.1 predicted. What stands in for a second method is
  replication of the same method — ten invocations of the bar and two of the whole cohort — which
  bounds sampling and run-to-run error but cannot detect a systematic error shared by every run.
  The tier-1 inversion is exactly the kind of systematic error a second predictor would be needed
  to confirm or refute.
- **The control calibrates dependence, not the ranking's validity.** That the 7OW6 receptor loses
  0.106 ipTM when the peptide is removed establishes the scale a design is judged against. It
  does not establish that a design matching that number binds the peptide the way the receptor
  does, and no design was shown to exceed it.
- **The wild-type pass was never replicated.** One invocation, so the tier-1 column has no
  between-run denominator — the single thing the peptide-dependence drop needed before its own
  verdicts could be trusted.
