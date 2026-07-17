# The geometry prior was measuring the wrong domain (2026-07-13)

touchstone's geometry tier z-scores a design's metal–donor bonds against an empirical prior. The
active prior was **CSD** — small-molecule metal–organic crystals. touchstone verifies **proteins**.

That mismatch is not academic. Scored against real, experimentally-determined metalloprotein sites,
the active prior **rejects reality**:

| metal | prior | real sites TRUSTED | real sites DEFERRED (off-manifold) |
|---|---|---|---|
| **Ni²⁺** | **CSD (was active)** | **41.3%** | **19.5%** |
| Ni²⁺ | MetalPDB (now default) | **96.6%** | 0.0% |
| **Co²⁺** | **CSD** | **41.4%** | 1.3% |
| Co²⁺ | MetalPDB | 96.6% | 0.3% |
| Cu²⁺ | CSD | 97.4% | 0.0% |
| Cu²⁺ | MetalPDB | 97.7% | 0.0% |

*(349 Ni / 309 Cu / 319 Co representative MetalPDB sites, first-shell N/O/S protein donors.)*

**One in five real nickel metalloprotein sites was called off-manifold**, and only 41% were trusted.
These are crystal structures of proteins that demonstrably bind nickel. A verifier that defers
reality isn't being conservative — it's miscalibrated.

### Cu²⁺ passes this table for the wrong reason

The Cu row above is the one to be careful with. The CSD Cu²⁺ prior is *not* fine — it is exactly as
small-molecule as the others. It scores 97.4% only because it is too blunt to reject anything:
Jahn–Teller elongation broadens it (std 0.211 Å) and its `cn_range` runs to **6**, while real protein
Cu²⁺ is **3–5 coordinate** (type-1 and type-2 copper are not octahedral).

Hand it an **octahedral CN6 copper site** — small-molecule/aqua geometry that no protein Cu site
adopts — and it returns `plausible (0.0σ, coordination in range)`. A *perfect* score. The protein
prior rejects it outright (`coordination number outside observed range`).

So the real-site table is a **sensitivity** check (does the prior accept real sites?), and a
wrong-domain prior can hide from it by being wide. What catches Cu is **specificity** (does it reject
wrong things?). This is the same blind spot that made the first version of the calibration test
useless — see the closing section — and it has teeth for the Cu work specifically: copper wants
square-planar / type-1 geometry, so **a Cu RLVR geometry reward built on the CSD prior would have
rewarded octahedral copper.** Pinned as `test_octahedral_copper_is_rejected`.

## Why

The CSD Ni²⁺ prior is centred at **2.074 Å**; real protein Ni–(N/O/S) bonds average **2.190 Å**. The
prior is **0.116 Å too short — more than 1σ of its own std (0.108 Å)**. Every protein-typical bond
therefore carries ~**+1.1σ of artificial strain** before any real strain is measured, against a
`trust_z` budget of 2.0σ. Over half the trust budget is spent on the domain error.

The corollary matters more than the miscalibration: the RLVR reward is built on this tier, so it was
**pulling BoltzGen toward 2.07 Å small-molecule chelate geometry when real nickel proteins want
2.19 Å.** The generator was being steered, precisely, at the wrong target.

## The fix, and what it costs

`metalpdb_reference.json` is now built and bundled ([`scripts/build_metalpdb_reference.py`](../../scripts/build_metalpdb_reference.py)),
and [`best_reference()`](../../src/touchstone/geometry/reference.py) already preferred MetalPDB over
CSD over PDB — so shipping the data flips the default. Nine metals (Ni/Cu/Co/Zn/Fe/Mn/Ca/Mg/Pt), up
from three. Two domain guards:

- **Solvent donors excluded — on both sides.** A designed structure carries no waters, so touchstone
  can only ever measure protein donors. Counting MetalPDB's waters shifts modal CN 6→4 for Ni/Co/Mn
  (they sit octahedrally with ~2 aqua ligands) and would mark a good 4-donor site with an open,
  water-fillable face as two donors short. The *parser* had the mirror-image hole —
  [`site_from_atoms`](../../src/touchstone/geometry/parse.py) filtered donors by element only, so a
  structure that does carry solvent (a crystal structure, a co-fold) had its water oxygens counted as
  donors against a prior that excludes them. Both sides now agree.
- **A 30-site floor.** Au³⁺ has 4 sites in the PDB and Pd²⁺ has 6. A mean±std from a handful of sites
  is noise with a decimal point, and a bundled prior is trusted silently by every verdict. Below the
  floor we ship nothing and report "no reference geometry" — absent beats fabricated. The builder also
  refuses to write an *empty* table: `best_reference()` selects on file existence, so an offline run
  would otherwise install a MetalPDB prior that `KeyError`s on every metal.

`MIN_CN = 3` drops the 58% of MetalPDB sites held by only 1–2 protein contacts (surface/adventitious
ions completed by water, not the buried sites a design targets). That selection moves the bond mean by
0.013 Å — immaterial against the 0.116 Å domain error — but it does mean `cn_range`'s lower bound is
that floor rather than a free observation.

## The same bug was in the precedent tier

Chasing the solvent question surfaced a second, independent instance.
[`metalpdb_precedents.json`](../../src/touchstone/data/metalpdb_precedents.json) — the motif→count
table behind the **default-on** [`PrecedentVerifier`](../../src/touchstone/geometry/precedent.py) — was
built with **no solvent filter at all**. Water-padded pseudo-motifs took over the counts:

| motif | with waters | protein donors only |
|---|---|---|
| `Ni-O6` (hexaaqua nickel — an ion sitting in solvent) | 22 | 7 |
| `Ni-N2O4` | 55 | 6 |
| **`Ni-S4`** (NiFe hydrogenase, Cys₄) | **3 → unprecedented** | **5 → trusted** |

With `min_hits=5`, touchstone was calling the **NiFe-hydrogenase nickel site unprecedented** — one of
the best-characterised nickel sites in biology — because aqua complexes had crowded out the real
motifs. Rebuilt with the shared filter, the Ni table now reads like protein chemistry (`Ni-N2O2`,
`Ni-N3`, `Ni-N4`) instead of aqua chemistry, and `Zn-S4` (zinc fingers) sits at 448 hits.

Both builders now import `SOLVENT_RESIDUES` from the parser rather than each keeping their own copy —
the two had already drifted (the reference builder's set was missing `H2O`).

**The honest cost: the geometry tier is now much more permissive.** On the RLVR pools:

| pool | CSD | MetalPDB |
|---|---|---|
| baseline (288) | 6.2% | 71.5% |
| round-1 (96) | 13.5% | 89.6% |
| round-2 (96) | 4.2% | 71.9% |
| round-3 (96) | 54.2% | 99.0% |

This is not a regression — it is the tier's true resolution. A 2σ gate on a correctly-centred prior
admits ~96% of real sites (textbook Gaussian: 95%), so it will also admit most geometrically-plausible
designs. **Real protein metal sites are genuinely heterogeneous** (Ni–N 2.173 ± 0.182 Å, Ni–O
2.239 ± 0.251 Å even split by donor element — the breadth is not Cys contamination, it is bidentate
carboxylates and backbone carbonyls). The CSD prior's apparent sharpness *was* its domain error:
it discriminated by rejecting protein geometry as such.

So geometry is a **weak filter**, and always was. Discrimination has to come from bond-valence, the
coordination tiers, and MLIP — which is exactly the defense-in-depth the consensus was designed
around. The prior that looked like the sharpest tier was borrowing its sharpness from the wrong domain.

## What survives

The RLVR arc replicates under the corrected prior — baseline 71.5% → round-1 89.6% → round-2 71.9%
(the MLIP-reward dip) → round-3 99.0%, the same shape as the documented CSD arc (5.9 → 21.9 → 6.3 →
66.7, [round-1 writeup](2026-07-01-rlvr-boltzgen-round1.md)). **The loop worked.** But the geometry
TRUST rates quoted throughout that writeup are CSD numbers, and the *magnitude* of the gains was
inflated by a yardstick that scored protein-like geometry as strained. Re-running the arc's reward
under the corrected prior is future work; the qualitative conclusion (the balanced reward lifts both
axes; the gains saturate) is unchanged.

## A follow-up the fix dissolved: Cu²⁺ "Jahn–Teller bimodality"

The [Cu prior audit](../specs/2026-07-12-cu-rlvr-selectivity.md) had flagged a second item: Cu²⁺'s
bond-length std is ~2× Ni's, read as the Jahn–Teller signature (a JT octahedron has ≈4 short + 2 long
bonds), with the recommendation to model it bimodally. **That was an artifact of the CSD prior.**

JT elongation is a *six-coordinate* effect, and in the protein domain only **3 of 309 Cu²⁺ sites
(1.0%)** are 6-coordinate — protein copper is 3–4 coordinate (type-1/type-2), so 99% of it
structurally cannot distort. The effect is real where it can occur (those 3 sites: 4 bonds at
2.107 Å, 2 at 2.630 Å — **+0.52 Å**, textbook) but it is not what makes the Cu prior wide. The width
is **donor heterogeneity** — Cu–N 2.119, Cu–O 2.230, Cu–S 2.273 Å — and the distribution is
broad-unimodal, not bimodal.

CSD *is* full of octahedral Cu chelates, so the JT reading was correct **there**. Carrying it into a
protein prior would have meant building a bimodal model of an effect 99% of the domain can't exhibit.
Worth noting on its own: fixing the domain didn't just correct the numbers, it **retired a planned
piece of work that only existed because of the wrong domain.**

## The pattern (fourth occurrence)

MetalHawk (confidently OOD on designed sites), MACE/OrbMol (confidently ranking metals they cannot
represent — [writeup](2026-07-13-mlip-cannot-rank-metals.md)), and now a prior confidently z-scoring a
domain it was never measured on. Every time, a component produced a plausible number and the pipeline
trusted it. The guard that works is the same one each time: **check the tier against ground truth it
must not fail** — MLIPs against Irving–Williams, geometry priors against real metalloprotein sites.

That check is now a test ([`TestShippedPriorCalibration`](../../tests/test_geometry.py)): 360 real
MetalPDB sites ([fixture](../../tests/fixtures/metalpdb_real_sites.json)), asserting the shipped prior
trusts ≥90% of them. Remove `metalpdb_reference.json` and it fails with *"trusts only 40% of real Ni2+
protein sites"*.

The first version of that test did **not** work: it probed a single site at the *mean* of the protein
distribution, where CSD happens to agree (1.1σ — strained but within the 2σ gate). It passed under the
broken prior. CSD's failure lives in the *spread* (std 0.108 vs the protein's 0.205) and in `cn_range`,
so a point probe at the centre is exactly where the bug hides. Caught in review — and worth recording,
because a green test asserting a calibration it never checked is the same failure as everything above,
one level up.
