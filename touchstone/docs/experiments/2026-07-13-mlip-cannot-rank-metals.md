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
[M(H₂O)₆]²⁺ and requires the **whole series** — stability rising Mn²⁺ < Fe²⁺ < Co²⁺ < Ni²⁺ < Cu²⁺,
then dropping back at Zn²⁺ — and the verifier `defer`s outright if it fails:

> `backbone cannot rank the divalent 3d series (fails Irving–Williams: stability must rise
> Mn<Fe<Co<Ni<Cu and drop at Zn) — it has no ligand-field/spin physics, so its metal ordering is
> meaningless`

**Requiring the series rather than the peak matters.** An `argmin == Cu2+` check is one bit, and one
bit is cheap to satisfy by accident: a backbone ranking Cu strongest but **Mn²⁺ second** — above Fe,
Co and Ni — would sail through it while still being physically meaningless. Pinned as
`CuPeakScrambledSpring` in [the tests](../../tests/test_selectivity.py).

With MACE and OrbMol both failing, **the selectivity tier is now inert** — which is the correct
behaviour. The `--selectivity` RLVR reward cannot steer a generator with a signal that ranks Mn²⁺
above Cu²⁺.

## Spin and charge follow the ion, not the element

Two silent-wrong-answer bugs sat under the gate, dormant only because no backbone passes it:

- **The d-count table was keyed by element**, so `Fe3+` (d⁵, 5 unpaired) was handed `Fe2+`'s spin
  (d⁶, 4). Likewise `Cu1+`→`Cu2+`, `Mn3+`→`Mn2+`. Now keyed by `(element, oxidation state)`.
- **The free-ion charge was hard-coded `+2`**, with the apo leg taking `total − 2`. For a trivalent
  ion both legs are then wrong — and Pd²⁺/Pt²⁺/**Au³⁺** are in the project's own metal panel, they
  are the e-waste recovery targets. `swap_metal` now moves element, oxidation state, charge and spin
  together.

And where the spin state genuinely isn't known — **Co³⁺ (d⁶) is high-spin in weak fields and low-spin
in most complexes** — the table says nothing and the tier defers, rather than substituting a singlet.
A wrong spin state doesn't crash; it just moves the ligand-field energy, which *is* the number this
tier reports. That is the same failure mode as everything else in this document, so it gets the same
answer: **refuse, don't guess.**

## The gap that remains: the probe is all-oxygen

The gate probes [M(H₂O)₆]²⁺ — hard, oxygen-only donors. The tier then judges **protein sites with N
and S donors**, and metal selectivity is largely a *donor-identity* (HSAB) effect: soft Cys-thiolate
and Met-thioether are exactly what a type-1 copper site uses to discriminate Cu from Ni/Co.

So a backbone that passes this gate has been shown to rank metals **on hard-donor sites only**. It
has not been shown to do so on the thiolate/imidazole chemistry the tier actually scores, and the
gate says nothing at all about metals outside the divalent 3d series (Pd²⁺/Pt²⁺/Au³⁺). Closing that
properly needs a soft-donor probe with reference data we don't have — so it is documented rather than
papered over. **Probe domain ≠ judged domain** is precisely the error that broke the geometry prior
([writeup](2026-07-13-geometry-prior-wrong-domain.md)); naming it here is cheaper than rediscovering it.

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
