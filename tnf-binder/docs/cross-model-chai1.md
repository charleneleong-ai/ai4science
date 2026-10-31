# Cross-model check: Chai-1 against Boltz-2

Every affinity and mouse-retention number in the submission came from Boltz-2, and BoltzGen
generated these designs using Boltz-2 weights — so that ranking was a model scoring its own
family's output. Chai-1 shares no weights with it. All 20 designs were re-folded against the
human and the mouse trimer, 5 diffusion samples each, 40 complexes.

Binder-to-target only, from `per_chain_pair_iptm` row X against chains A/B/C. The global number
is unusable for the same reason it was unusable for Boltz-2: the trimer's own A–B/B–C/A–C
interfaces score **0.80–0.81** and swamp the binder's, pulling `interface_ptm` to 0.719 on a
complex whose binder interface is 0.11. **That artefact now reproduces independently in a second
model**, which retroactively justifies the correction rather than merely asserting it.

## The retention ordering has no independent support

| quantity | Spearman, Chai-1 vs Boltz-2 |
| --- | --- |
| **mouse retention** | **+0.046** |
| human binder-to-target confidence | +0.344 |

A rank correlation of +0.046 is no relationship at all. Mouse retention is the term that
separates designs *within* each switch tier in the submission ranking, and it is not reproduced
by an independent model.

The clearest disagreement is the free_1 pair, the two designs most argued over:

| design | Boltz-2 retention | Chai-1 retention |
| --- | --- | --- |
| tnf2_09 | 0.466 — **worst of 20** | 1.057 — **best of 20** |
| tnf2_10 | 0.482 — second worst | 1.056 — second best |

A complete inversion on exactly the rows whose placement was debated most.

## Boltz-2's discrimination is not reproduced

| model | mean | range |
| --- | --- | --- |
| Chai-1 | 0.117 | 0.106 – 0.159 |
| Boltz-2 | 0.207 | 0.151 – 0.434 |

Chai-1 places all 20 designs in a 0.05-wide band and sees essentially no difference between
them. Boltz-2's apparent spread — especially the free_1 pair at 0.424 and 0.434, which was read
as "roughly twice the affinity of the helix_0 series" — is not reproduced. Under Chai-1 the
free_1 pair is unremarkable, at 0.111 and 0.113 against a pack mean of 0.117.

## What this changes

The switch-tier ordering (tier A above tier B above tier C) rests on pKa evidence verified three
ways and stands. The ordering *within* each tier rests on mouse retention, and that term now has
no independent support. It should be reported as a Boltz-2 reading rather than a measurement,
and the claim that the free_1 pair carries materially higher affinity should be withdrawn: one
model says so, another says it is ordinary, and nothing experimental has been measured.

This landed after the Challenge 2 deadline. It is recorded here for the result and because the
same check is a standing requirement for subsequent challenges.
