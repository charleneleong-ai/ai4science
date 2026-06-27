"""Run a structure through the verifier stack → one JSON-able verdict.

The shared engine behind the CLI and the MCP server. Runs every verifier that can
judge a structure inline, with graceful degradation:

  - **always** (pure-Python, instant): geometry z-score + bond-valence.
  - `deep=True`: MLIP relaxation + MLIP-MD (need a GPU backend; skipped if absent).

Stages that need an external input a bare structure can't supply — Mogul (a CSD
licence), co-fold (a second predictor's structure), expression (a sequence scorer),
global thermostability (an MD/Tm scorer) — are reported under `not_run` rather than
guessed. A verifier that can't run (no backend) is *skipped* (excluded from the
verdict), distinct from one that ran and *deferred*.

Consensus is defense-in-depth: `trust` only if every verifier that ran trusts, `defer`
if any defers (or none ran), else `weak`.
"""

from __future__ import annotations

from pathlib import Path

from .core import BinderDesign, Verdict, element_symbol
from .geometry.bond_valence import BondValenceVerifier
from .geometry.parse import coordination_site
from .geometry.reference import PDBReference
from .geometry.verifier import GeometryVerifier
from .physics.mlip import MLIPDynamicsVerifier, MLIPVerifier, make_backbone  # light import; heavy load lazy

# stages with a library verifier but no inline-available input (licence / prediction /
# scorer) — advertised so an agent knows the full stack and how to enable each.
_NEEDS_INPUT = {
    "mogul": "a CSD licence (Mogul / CSD Python API)",
    "cofold": "a co-fold prediction (scripts/chai_crosscheck or allmetal3d_crosscheck)",
    "expression": "a sequence scorer (scripts/expression_score)",
    "thermostability": "an MD/Tm scorer (scripts/thermostability_score)",
}


def _as_dict(v: Verdict) -> dict:
    return {"label": v.label, "score": round(v.score, 3), "trust": v.trust, "ood": v.ood, "reason": v.reason}


def verify_structure(structure: str | Path, metal: str = "Ni2+", deep: bool = False, cutoff: float = 2.8) -> dict:
    """Verify a metal-coordination structure. Returns per-verifier verdicts, a `not_run`
    map of stages needing inputs, and a trust/weak/defer consensus."""
    site = coordination_site(structure, element_symbol(metal).upper(), metal, cutoff)
    design = BinderDesign("", site, generator="external", generator_confidence=0.0, source=str(structure))

    verifiers = {"geometry": GeometryVerifier(PDBReference()), "bond_valence": BondValenceVerifier()}
    results: dict[str, dict] = {}
    if deep:
        try:  # build the backbone once; share it across both MLIP verifiers
            calc = make_backbone("mace_mp")
        except Exception as e:  # no GPU backend ⇒ skip both (not a defer that tanks consensus)
            for n in ("mlip", "mlip_md"):
                results[n] = {"skipped": f"no MLIP backend: {type(e).__name__}"}
        else:
            verifiers["mlip"] = MLIPVerifier(calculator=calc)
            verifiers["mlip_md"] = MLIPDynamicsVerifier(calculator=calc)

    counted: list[str] = []
    for name, verifier in verifiers.items():
        try:
            verdict = verifier.verify(design)
        except Exception as e:  # unexpected per-verifier failure ⇒ skipped, not counted
            results[name] = {"skipped": f"{type(e).__name__}: {e}"}
            continue
        results[name] = _as_dict(verdict)
        counted.append(verdict.label)

    consensus = (
        "defer" if (not counted or "defer" in counted)
        else "trust" if all(label == "trust" for label in counted)
        else "weak"
    )
    return {
        "structure": str(structure),
        "metal": metal,
        "coordination_number": site.coordination_number,
        "donors": list(site.ligand_elems),
        "verifiers": results,
        "not_run": _NEEDS_INPUT,
        "consensus": consensus,
    }
