---
description: Verify a designed metal-binding structure with touchstone and report the trust/weak/defer verdict.
argument-hint: <structure.pdb|cif> [metal — defaults to Ni2+]
allowed-tools: mcp__touchstone__verify_metal_binder, Bash(touchstone:*), Read
---

Verify the metal-binding design at `$1` against the touchstone verifier stack. Use the metal label in `$2` if given, otherwise `Ni2+`.

1. Call `verify_metal_binder(structure_path="$1", metal=<$2 or "Ni2+">)` — or `touchstone verify "$1" --metal <metal>`.
2. Report the **consensus** (trust/weak/defer), the per-verifier verdicts, and whether it clears the wet-lab bar (only `trust` does).
3. If `weak`/`defer`, name which verifier flagged it (geometry σ, bond-valence Δ, co-fold disagreement, …) and what to change before re-verifying.
