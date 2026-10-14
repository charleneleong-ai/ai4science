# tnf-binder — a conditional TNF-α binder, and the verifiers that gate it

Adaptyv × Anthropic Protein Design Competition, challenge 02. One binder, three objectives:
bind human TNF-α, cross-react with mouse TNF-α, and bind at pH 7.4 while releasing at pH 6.0.

This directory is a worked instance of the repo's through-line — **the generator is the
commoditised part; the verifier is where reliability lives.** Four generators were
interchangeable here. The verifiers decided everything.

## The headline: what the verifiers caught

| verifier | what it caught |
| --- | --- |
| sequence novelty (exact k-mer vs PDB + Swiss-Prot) | the **top-ranked construct was verbatim ubiquitin** — best interface score, lowest interface PDE, an ideal-pKa histidine, and disqualified on sight |
| sequence composition | the top of a 40-design shortlist was **poly-alanine** (A=59–76%); 54 of 57 constructs failed |
| cross-model agreement | AF2 scored its own optimised designs at **i_pTM 0.92**; an independent model scored the same structures **0.21** |
| fold novelty (FoldSeek + TM-align) | predicted the one **real-world novelty failure** correctly — a β-sandwich matching CATH 2.60.40 at TM 0.893 |
| PROPKA pKa prediction | the engineered histidines from first-principles reasoning were **all dead** (pKa 3.1–5.7 against a 6.0–7.4 window) |
| [sampling depth](docs/sampling-depth.md) (20 diffusion samples, not 3) | **two of the four switch positions have no determined pKa** — His21 and His35 span 4.3 and 5.7 pH units across structures the folding model scores equally, and three samples hid it |

Every one of those passed the checks that preceded it. None would have been caught by a
generator, and none by eyeballing a confidence score.

## The empirical calibration (n=15)

Proteinbase's automatic novelty check is the only external measurement this campaign gets. Of
15 designs scored across two submissions, **one failed** and 14 passed: submission 1 put up 3
(one failure, two at 3/4), submission 2 put up 12 and **all 12 passed** — eleven at 3/4 and
`free1_T73H` at **4/4**.

The fold gate predicts the one failure correctly and **under-scores every pass**, because its
search measures more similarity than Proteinbase's does. Treat its verdict as a floor, never as
the score:

| design | real outcome | fold verdict | TM | coverage | predicted |
| --- | --- | --- | --- | --- | --- |
| `V__pocket_free_3__A73H` | **failed** | high | 0.893 | 0.888 | 2 — fail ✓ |
| `V__pocket_helix_0__A35H` | passed 3/4 | moderate | 0.776 | 0.802 | 2 — under-scored |
| `P__pocket_helix_0` | passed 3/4 | moderate | 0.759 | 0.802 | 2 — under-scored |
| `free1_T73H` | passed **4/4** | moderate | 0.704 | 0.973 | 3 — under-scored |

An earlier claim that this reproduced all three submission-1 verdicts was **wrong**: the two
passes scored correctly only because a missing sequence-identity value defaulted to zero. The
tests now encode the disagreement rather than hiding it.

The failing design's closest match is CATH **2.60.40**, the immunoglobulin-like β-sandwich.
Sequence identity was near-identical across all three (37.5% vs 38.1%), so **fold similarity —
not sequence novelty — was the binding constraint.** An exact-k-mer screen cannot see this.

**Novelty and the objectives pulled in opposite directions.** The campaign's only 4/4,
`free1_T73H`, is also its worst design on the two things ranked higher: 47% of its human
interface retained against mouse, and 0 of its 1 histidine holding the pKa window across
diffusion samples. Novelty is a gate to pass (≥3/4), not a score to maximise.

**One residue is worth one novelty point, and our gate cannot see it.** `free1_T73H` and
`free1_parent` are the same 74-residue backbone one substitution apart, with the same nearest
hit (O27230, 35.8% identity) and the same fold verdict. Proteinbase scored them **4/4 and
3/4**; our sequence-only gate scored both 3. So the gate is not merely conservative — it has
no resolving power at the 3/4–4/4 boundary, and whatever moved is invisible to identity and
TM-score alike. With no model of that component and one observation, chasing 4/4 would be a
blind search over single substitutions, paid for one upload at a time, for no gain on any
ranked objective.

## Submission 2: the two measurements that reordered everything

Both came from re-measuring something already treated as known, and both overturned a ranking
that had already been built — the per-design numbers are in
[`docs/submission-v2-evidence.csv`](docs/submission-v2-evidence.csv).

**Every pKa came from a single predicted structure.** Boltz writes three diffusion samples per
design, so each histidine can be scored three times against three independent guesses at the
same structure. Re-scored that way, **8 of 23 histidines leave the window** — `R21H` alone
ranges pKa 3.99–6.70 where its design structure said 7.05 — and only five designs hold every
histidine in window on all three samples. [`switch.SampledHis`](src/tnf_binder/switch.py) is
that measurement; the ranking statistic is `SwitchProfile.delta_q_floor`, the charge the
interface gains between pH 7.4 and 6.0 with each histidine taken at its *least favourable*
sample. Distance from pKa 6.7 is what costs — a pKa of 7.6 is as bad a bet as one of 5.8.

**Mouse cross-reactivity was being inferred from a proxy, and the proxy was wrong.** Folding
all 12 designs against the mouse trimer (UniProt **P06804**) leaves 9/12 retaining ≥70% of their
human interface, median 86% — but the `free_1` family collapses to **47%**, despite carrying the
campaign's *highest* conserved-epitope coverage (0.79). The design promoted on that proxy was
the one the proxy was wrong about.

## The science that held up

The epitope was computed, not taken from literature: TNF residues within 4.5 Å of TNFR2 in
**3ALQ**, intersected with a human/mouse alignment (124/157 = 79.0% identity, reproducing the
stated figure). That gives 33 contact positions, 24 conserved and 9 divergent.

**R31→Q is the trap.** It carries the largest SASA in the epitope (193 Å²) and sits beside R32
in the 29–33 loop that adalimumab and infliximab engage — which is why neither drug is mouse
cross-reactive. Any design gripping that loop wins objective 1 and loses objective 2.

## What is here

| path | what |
| --- | --- |
| [`src/tnf_binder/`](src/tnf_binder) | the package: the target and its numbering, the four verifiers, confidence reading, and the `tnf-binder` CLI (`epitope`, `verify`, `validate`, `confidence`, `wandb-log`) |
| [`scripts/`](scripts) | the design-box pipeline, where the heavy dependencies live. Ten entry points, each a typer app: [`gen.py`](scripts/gen.py) (`spec`/`screen`) writes BoltzGen specs, [`his.py`](scripts/his.py) (`scan`/`combos`/`campaign`) runs the PROPKA+PyRosetta histidine work, [`metrics.py`](scripts/metrics.py) (`collect`/`join`) reads fold confidence, [`submit.py`](scripts/submit.py) (`shortlist`/`rank`/`final`) gates and assembles an upload |
| [`specs/`](specs) | every generator input that ran, one directory per tool (the JSON/YAML split is the tools', not ours): [`bindcraft1/`](specs/bindcraft1) the abandoned first campaign, [`bindcraft2/`](specs/bindcraft2) three dual-ortholog rounds that scored zero, [`boltzgen/`](specs/boltzgen) the screens that produced the submission plus a 2x2 over method and target |
| [`docs/method.md`](docs/method.md) | epitope choice, the pH mechanism, scoring definitions |
| [`docs/experiment-card.md`](docs/experiment-card.md) | the stack, the funnel with denominators, and what is *not* established |
| [`submission/`](submission) | the ranked CSV as submitted |

## Running it

```bash
cd tnf-binder
uv run --locked --extra dev pytest -q          # 202 tests, no GPU and no PyRosetta needed
uv run tnf-binder epitope                      # the computed epitope and its divergent positions
uv run tnf-binder validate submission/tnf_submission_v2.csv
```

The verifiers, the pKa arithmetic and the submission gates are all pure and tested here. The
generators (BoltzGen, BindCraft2), the folding (Boltz-2) and PROPKA/PyRosetta run on a GPU box
and are driven by [`scripts/`](scripts); nothing in the package imports them.

## Honest status

*Status as of 2026-10-11.* Two submissions scored, 15 designs, 14 novelty passes, one 4/4.
No design here is a validated binder. Interface confidence is **model-dependent and
unresolved**: Boltz-2 and BoltzGen's folding step disagree by 2–3× on identical structures, and
AF2 disagrees with both. A generated protein is a candidate, not a result — see
[`docs/experiment-card.md`](docs/experiment-card.md) for what would be needed to make it one.
