# Can touchstone tell Ni from Cu/Co? Not from geometry alone — and that's the point

**2026-06-25** · Selectivity — binding Ni²⁺ *over* the Cu²⁺/Co²⁺/Fe²⁺ also in leachate — is the
real critical-mineral-recovery problem. The honest result: **the geometry verifier cannot deliver
it**, and saying so clearly is more useful than a fabricated selectivity number.

## Geometry can't discriminate

The three metals' empirical PDB geometries overlap heavily:

| metal | bond mean ± std (Å) | CN range |
| --- | --- | --- |
| Ni²⁺ | 2.154 ± 0.183 | 4–6 |
| Cu²⁺ | 2.144 ± 0.179 | 3–5 |
| Co²⁺ | 2.165 ± 0.186 | 3–6 |

`selectivity_profile(design, verifier, ["Ni2+", "Cu2+", "Co2+"])` re-scores a site as each metal.
Across the Ni-trusted designs: **0/8 were Ni-selective** — every one also `trust`s as Cu and Co,
because a CN4 site at ~2.15 Å sits in-range for all three. Geometry only discriminates at the CN
*extremes* (a CN6 site excludes Cu, whose range tops at 5; a CN3 site excludes Ni). In the CN4–5
band where real designs land, it can't tell them apart.

## Physics doesn't rescue it cheaply either

The obvious next move — swap the metal in a trusted site and let xtb say which it prefers — is
**artifact-prone**, not a clean answer:

- an *isolated* coordination cluster (no protein scaffold) over-tightens on optimisation (bonds
  collapse to ~1.8 Å) and loses donors, so the relaxed geometry isn't the in-protein geometry;
- GFN2-xTB is unreliable for transition-metal **spin states / d-electron energetics** (Cu²⁺ d⁹,
  Co²⁺ d⁷ are open-shell), so the binding *energies* aren't trustworthy.

A defensible selectivity calculation needs constrained full-system simulation or higher-level QM —
not a quick cluster optimisation.

## The lesson

Selectivity has to be **designed in**, not verified after the fact:

- **donor identity (HSAB)** — soft S-donors (Cys) discriminate soft/borderline metals; geometry
  ignores donor type entirely;
- **rigid, pre-organised cavities** sized for one metal's ionic radius and coordination preference.

The verifier's honest role here is twofold: provide the multi-metal profile (`selectivity_profile`),
and **flag that "geometry-plausible ≠ selective"** so a generator isn't trusted to have solved
selectivity when it has only produced a plausible divalent-metal site. Knowing the limit of the
geometry oracle is itself a result.
