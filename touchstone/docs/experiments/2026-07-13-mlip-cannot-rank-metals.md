# MLIPs cannot rank the divalent 3d series — the selectivity tier must refuse (2026-07-13)

The MLIP metal-swap selectivity tier ([`MLIPSelectivityVerifier`](../../src/touchstone/physics/selectivity.py))
asks *"does this site prefer its target metal over the competitors?"* — the discrimination geometry
can't make, and the lever the Cu RLVR reward was built on
([spec](../specs/2026-07-12-cu-rlvr-selectivity.md)).

**It doesn't work, and it can't be made to work with any MLIP backbone we have.** This documents
how we found out, why, and the fix — and retracts the findings that were built on it.

## How it surfaced

Scoring the 96 round-4 Ni designs with the selectivity tier gave a striking result: **0/96
Ni-selective**, with **Co²⁺ preferred in 89/96** (median margin −0.59 eV). That looked like a major
finding — the RLVR rounds had optimised *binding* and quietly produced designs that prefer a
different metal.

It was too clean. A site's metal preference follows the **Irving–Williams series** —
Mn < Fe < Co < **Ni < Cu** > Zn, with **Cu²⁺ the peak** — which holds for essentially any N/O donor
set and is among the most robust empirical trends in coordination chemistry. Co²⁺ beating both Ni²⁺
and Cu²⁺ across 192 structurally diverse designs is not a property of the designs. It is a property
of the potential.

## The test: [M(H₂O)₆]²⁺, the system the series was measured on

| backbone | ranking (strongest → weakest) | Cu²⁺ the peak? |
|---|---|---|
| **Irving–Williams (truth)** | Cu > Ni > Co > Fe > Mn (Zn below Cu) | — |
| **MACE-MP** | Mn > Co > Fe > **Cu** > Zn > Ni | ❌ (Cu 4th, Ni last) |
| **OrbMol** | Fe > Mn > Zn > Co > Ni > **Cu** | ❌ (Cu **last**) |

Both were given the fair version of the test: correct **+2 charge**, correct **per-metal spin**
(high-spin d-count: Mn²⁺ d⁵ → 5 unpaired … Zn²⁺ d¹⁰ → 0), and a **full relaxation** so each metal
finds its own M–O distance (a rigid probe would penalise Cu²⁺, which is Jahn–Teller distorted).
They still fail.

## Why — and why it is not fixable by tuning

Metal preference is set by **ligand-field stabilisation**: d-electron count and spin state. Two
independent failures compound:

1. **MACE-MP carries no charge or spin state at all.** It is trained on inorganic crystals. Handed a
   bare Ni/Mn/Cu atom in a cluster it returns a number, but that number cannot encode the physics
   that decides the ordering — it appears to track ionic size / cohesive energy instead. It is
   *confidently wrong*, in exactly the way [MetalHawk](2026-07-07-metalhawk-ood-designed-sites.md)
   was.
2. **The tier's own swap was incoherent.** `swap_metal` changed the element symbol and nothing else,
   so every metal in the panel was evaluated at the **same spin multiplicity** — comparing Mn²⁺ (d⁵)
   and Zn²⁺ (d¹⁰) as if they were the same electronic system. Fixed here (spin now follows the
   metal), but fixing it does not rescue either backbone.

OrbMol is the interesting case: it is a *molecular* potential and **refuses to run without charge and
spin** (`ValueError: atoms.info must contain both 'charge' and 'spin'`) — the model that knows this
physics is required declines to guess. It still cannot rank the series.

## The fix: a validity gate, not a better number

A tier that cannot reproduce a known reference trend has no business emitting verdicts.
[`ranks_irving_williams`](../../src/touchstone/physics/selectivity.py) probes the backbone on
[M(H₂O)₆]²⁺ and requires Cu²⁺ at the peak; the verifier `defer`s outright if it fails:

> `backbone cannot rank the divalent 3d series (fails Irving–Williams: Cu2+ is not the peak) — it
> has no ligand-field/spin physics, so its metal ordering is meaningless`

With MACE and OrbMol both failing, **the selectivity tier is now inert** — which is the correct
behaviour. The `--selectivity` RLVR reward cannot steer a generator with a signal that ranks Mn²⁺
above Cu²⁺.

## Retractions

Everything derived from the ungated tier is void:

- ❌ "0/96 of the RLVR Ni designs are Ni-selective"
- ❌ "Co²⁺ out-competes Ni²⁺ in 89/96 designs (median −0.59 eV)"
- ❌ the Ni-motif vs Cu-motif (hard- vs soft-donor) head-to-head

The Cu-motif experiment is **untested**, not failed — the yardstick was broken.

## What would actually work

Metal selectivity needs a method that carries ligand-field physics: **DFT with explicit spin states**
(the reference answer, slow), or a semi-empirical method with d-orbital treatment. It is not an MLIP
tier. Until a backbone passes the gate, touchstone reports selectivity as `defer` rather than
guessing — and the gate is the mechanism that keeps that honest.

## The pattern (third occurrence)

A fast surrogate, used outside its training distribution, produces a plausible number that a
pipeline then trusts. MetalHawk (confidently OOD on designed sites), and now MACE/OrbMol
(confidently ranking metals they cannot represent). The lesson generalises beyond this tier:
**the verifier needs validity gates on its own tiers, not only on the designs it judges.**
