"""Run a structure through the verifier stack → one JSON-able verdict.

The shared engine behind the CLI and the MCP server. Runs every verifier that can
judge a structure inline, with graceful degradation:

  - **always** (pure-Python, instant): geometry z-score + bond-valence.
  - `deep=True`: MLIP relaxation + MLIP-MD (need a GPU backend; skipped if absent).
    Protonates the structure first (OpenBabel) so MACE sees a chemically complete
    site — skipped, with no protonation, if OpenBabel isn't installed.

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

from .cofold import CofoldCrossCheck
from .core import BinderDesign, Verdict, element_symbol
from .geometry.bond_valence import BondValenceVerifier
from .geometry.coordination import CoordinationGeometryVerifier, CoordinationSymmetryVerifier
from .geometry.metalhawk import MetalHawkVerifier
from .geometry.parse import coordination_site
from .geometry.reference import best_reference
from .geometry.verifier import GeometryVerifier
from .physics.mlip import MLIPDynamicsVerifier, MLIPVerifier, make_backbone  # light import; heavy load lazy
from .pipeline import stress_profile

# the lightweight verifiers are stateless over read-only package data — build once and reuse,
# so a batch `rank` doesn't re-parse JSON per design. Geometry uses the sharpest reference on
# hand: the CSD metal–organic prior if it's been built, else the committed PDB reference.
_REFERENCE = best_reference()
_GEOMETRY = GeometryVerifier(_REFERENCE)
_BOND_VALENCE = BondValenceVerifier()
_COORD_SYMMETRY = CoordinationSymmetryVerifier()  # nVECSUM: is the metal enclosed?
_COORD_GEOMETRY = CoordinationGeometryVerifier()  # polyhedron shape vs ideal

# stages with a library verifier but no inline-available input (licence / prediction /
# scorer) — advertised so an agent knows the full stack and how to enable each.
_NEEDS_INPUT = {
    "mogul": "a CSD licence (Mogul / CSD Python API)",
    "metalhawk": "MetalHawk geometry predictions (scripts/metalhawk_score.py — open, no licence)",
    "trs": "an apo (unbound) structure to diff against (topology-reorganization on binding)",
    "cofold": "a co-fold prediction (scripts/chai_crosscheck, alphafold3_crosscheck, or allmetal3d_crosscheck)",
    "expression": "a sequence scorer (scripts/expression_score)",
    "thermostability": "an MD/Tm scorer (scripts/thermostability_score)",
}


def _as_dict(v: Verdict) -> dict:
    d = {"label": v.label, "score": round(v.score, 3), "trust": v.trust, "ood": v.ood, "reason": v.reason}
    if v.metrics:  # machine-readable numbers behind the verdict (σ, BVS, drift…) when the verifier exposes them
        d["metrics"] = v.metrics
    return d


# the full verifier stack in cost order — the unified `stack` view lists every tier with its
# status so the consensus is auditable, even tiers that didn't run on this input
_STACK_ORDER = (
    "geometry", "bond_valence", "coord_symmetry", "coord_geometry", "metalhawk",
    "mogul", "mlip", "mlip_md", "trs", "cofold", "expression", "thermostability",
)


def _stack(results: dict) -> list[dict]:
    """One entry per stack tier, in cost order: ran (with its verdict), skipped (backend
    absent), or needs_input (licence / scorer / co-fold) — the complete per-stage picture."""
    rows = []
    for stage in _STACK_ORDER:
        if stage in results and "label" in results[stage]:
            rows.append({"stage": stage, "status": "ran", **results[stage]})
        elif stage in results:  # ran-but-skipped (e.g. MLIP with no backend)
            rows.append({"stage": stage, "status": "skipped", "detail": results[stage]["skipped"]})
        elif stage in _NEEDS_INPUT:
            rows.append({"stage": stage, "status": "needs_input", "detail": _NEEDS_INPUT[stage]})
        elif stage in ("mlip", "mlip_md"):  # only attempted with deep=True + a GPU backend
            rows.append({"stage": stage, "status": "needs_input", "detail": "pass deep=True (needs a GPU backend)"})
    return rows


_AUTO = object()  # verify_structure default: build the MLIP backbone per call (a batch passes a shared one)


def mlip_backbone():
    """The default MLIP backbone (MACE-MP), or None if no backend is installed — so the
    caller can skip the MLIP tier cleanly rather than deferring. Build it once and pass it
    to a batch of verify_structure calls to avoid reloading the model per design."""
    try:
        return make_backbone("mace_mp")
    except Exception:  # no GPU/torch backend ⇒ MLIP tier unavailable
        return None


def verify_structure(
    structure: str | Path, metal: str = "Ni2+", deep: bool = False, cutoff: float = 2.8, stress: bool = False,
    cofold_provider=None, metalhawk_scorer=None, calc=_AUTO,
) -> dict:
    """Verify a metal-coordination structure. Returns per-verifier verdicts, a `not_run`
    map of stages needing inputs, and a trust/weak/defer consensus. With `stress`, also
    re-verify the site under extreme-condition perturbations (acidic-leachate bond stretch,
    low-pH donor protonation) → a `stress` map {neutral/leachate/low_pH: verdict}: does it
    hold up in the real recovery process? `cofold_provider` (a design → predicted
    CoordinationSite callback over Chai-1 / AllMetal3D outputs) adds the independent co-fold
    cross-check tier; `metalhawk_scorer` (a design → MetalHawkPrediction over
    scripts/metalhawk_score.py output) adds the open geometry-distortion tier. `calc` is an
    internal knob for batch callers (`rank_structures`) to share one MLIP backbone; leave it
    default."""
    site = coordination_site(structure, element_symbol(metal).upper(), metal, cutoff)
    design = BinderDesign("", site, generator="external", generator_confidence=0.0, source=str(structure))

    verifiers = {
        "geometry": _GEOMETRY,
        "bond_valence": _BOND_VALENCE,
        "coord_symmetry": _COORD_SYMMETRY,
        "coord_geometry": _COORD_GEOMETRY,
    }
    if cofold_provider is not None:  # independent predictor (Chai-1 / AllMetal3D) corroboration
        verifiers["cofold"] = CofoldCrossCheck(cofold_provider)
    if metalhawk_scorer is not None:  # independent ANN geometry-distortion oracle
        verifiers["metalhawk"] = MetalHawkVerifier(metalhawk_scorer)
    results: dict[str, dict] = {}
    if deep:
        if calc is _AUTO:  # single call ⇒ build per call; a batch hands in a shared backbone (or None)
            calc = mlip_backbone()
        if calc is None:  # no backend ⇒ skip both (not a defer that tanks consensus)
            for n in ("mlip", "mlip_md"):
                results[n] = {"skipped": "no MLIP backend (install touchstone[mace])"}
        else:  # share the one backbone across both MLIP verifiers (they protonate internally)
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
    result = {
        "structure": str(structure),
        "metal": metal,
        "coordination_number": site.coordination_number,
        "donors": list(site.ligand_elems),
        "reference": _REFERENCE.source,  # which geometry prior backed the z-score (CSD or PDB)
        "verifiers": results,
        "not_run": _NEEDS_INPUT,
        "stack": _stack(results),  # full per-tier breakdown (ran / skipped / needs_input), cost order
        "consensus": consensus,
    }
    if stress:  # robustness map: does the site hold its verdict across the operating envelope?
        # geometry-tier by design (independent of `deep`): the perturbations are geometric
        # (bond stretch, donor protonation), so the z-score is the natural judge — and running
        # the MLIP tier across every condition would multiply GPU cost for little extra signal.
        result["stress"] = {cond: _as_dict(v) for cond, v in stress_profile(design, _GEOMETRY).items()}
    return result
