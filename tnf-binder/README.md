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

Every one of those passed the checks that preceded it. None would have been caught by a
generator, and none by eyeballing a confidence score.

## The one empirical calibration

One design was submitted and **failed** Proteinbase's automatic novelty check; two others
passed. The fold gate predicts the failure correctly, and is **conservative** on the passes —
it under-scores them, because its search measures more similarity than Proteinbase's does:

| design | real outcome | fold verdict | TM | coverage | predicted score |
| --- | --- | --- | --- | --- | --- |
| `V__pocket_free_3__A73H` | **failed** | high | 0.893 | 0.888 | 2 — fail ✓ |
| `V__pocket_helix_0__A35H` | passed | moderate | 0.776 | 0.802 | 2 — under-scored |
| `P__pocket_helix_0` | passed | moderate | 0.759 | 0.802 | 2 — under-scored |

An earlier claim that this reproduced all three verdicts was **wrong**: the two passes scored
correctly only because a missing sequence-identity value defaulted to zero. The tests now
encode the disagreement rather than hiding it.

The failing design's closest match is CATH **2.60.40**, the immunoglobulin-like β-sandwich.
Sequence identity was near-identical across all three (37.5% vs 38.1%), so **fold similarity —
not sequence novelty — was the binding constraint.** An exact-k-mer screen cannot see this.

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
| [`scripts/`](scripts) | the design-box pipeline, where the heavy dependencies live: spec generation (BoltzGen), the PROPKA/PyRosetta His scans, batch scoring, submission assembly |
| [`specs/`](specs) | BoltzGen design specs, BindCraft/BindCraft2 campaign settings |
| [`docs/method.md`](docs/method.md) | epitope choice, the pH mechanism, scoring definitions |
| [`docs/experiment-card.md`](docs/experiment-card.md) | the stack, the funnel with denominators, and what is *not* established |
| [`submission/`](submission) | the ranked CSV as submitted |

## Honest status

No design here is a validated binder. Interface confidence is **model-dependent and
unresolved**: Boltz-2 and BoltzGen's folding step disagree by 2–3× on identical structures, and
AF2 disagrees with both. A generated protein is a candidate, not a result — see
[`docs/experiment-card.md`](docs/experiment-card.md) for what would be needed to make it one.
