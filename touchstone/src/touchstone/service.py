"""Run a structure through the verifier stack → one JSON-able verdict.

The shared engine behind the CLI and the MCP server. **Lightweight by default**
(geometry z-score + bond-valence — instant, no GPU); `deep=True` adds the MLIP
relaxation when the optional backend is installed. Co-fold and MD verifiers run via
the dedicated inference scripts (they need a GPU box / precomputed predictions), so
they are out of a single inline call's scope.

Consensus is defense-in-depth: `trust` only if every verifier trusts, `defer` if any
defers (or none could run), else `weak`.
"""

from __future__ import annotations

from pathlib import Path

from .core import BinderDesign, Verdict, element_symbol
from .geometry.bond_valence import BondValenceVerifier
from .geometry.parse import coordination_site
from .geometry.reference import PDBReference
from .geometry.verifier import GeometryVerifier
from .physics.mlip import MLIPVerifier  # light at import-time; heavy ase/mace load lazily on use


def _as_dict(v: Verdict) -> dict:
    return {"label": v.label, "score": round(v.score, 3), "trust": v.trust, "ood": v.ood, "reason": v.reason}


def verify_structure(structure: str | Path, metal: str = "Ni2+", deep: bool = False, cutoff: float = 2.8) -> dict:
    """Verify a metal-coordination structure. Returns per-verifier verdicts + a
    trust/weak/defer consensus, ready to serialise for an agent or the CLI."""
    site = coordination_site(structure, element_symbol(metal).upper(), metal, cutoff)
    design = BinderDesign("", site, generator="external", generator_confidence=0.0, source=str(structure))

    verifiers = {"geometry": GeometryVerifier(PDBReference()), "bond_valence": BondValenceVerifier()}
    if deep:
        verifiers["mlip"] = MLIPVerifier(backbone="mace_mp")  # lazy; errors if no GPU backend

    results: dict[str, dict] = {}
    for name, verifier in verifiers.items():
        try:
            results[name] = _as_dict(verifier.verify(design))
        except Exception as e:  # an unavailable/failed verifier is reported, not fatal
            results[name] = {"error": f"{type(e).__name__}: {e}"}

    labels = [r["label"] for r in results.values() if "label" in r]
    consensus = (
        "defer" if (not labels or "defer" in labels)
        else "trust" if all(label == "trust" for label in labels)
        else "weak"
    )
    return {
        "structure": str(structure),
        "metal": metal,
        "coordination_number": site.coordination_number,
        "donors": list(site.ligand_elems),
        "verifiers": results,
        "consensus": consensus,
    }
