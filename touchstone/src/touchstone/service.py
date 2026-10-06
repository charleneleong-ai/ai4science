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

from collections.abc import Callable
from pathlib import Path
from typing import Literal

from .cofold import CofoldCrossCheck
from .core import BinderDesign, Verdict, Verifier, element_symbol
from .expression import ExpressionVerifier
from .geometry.bond_valence import BondValenceVerifier
from .geometry.coordination import CoordinationGeometryVerifier, CoordinationSymmetryVerifier
from .geometry.metalhawk import MetalHawkVerifier
from .geometry.mogul import MogulVerifier
from .geometry.parse import coordination_site
from .geometry.precedent import MotifSelectivityVerifier, PrecedentVerifier
from .geometry.reference import best_reference
from .geometry.verifier import GeometryVerifier
from .physics.mlip import MLIPDynamicsVerifier, MLIPVerifier, make_backbone  # light import; heavy load lazy
from .physics.selectivity import MLIPSelectivityVerifier
from .physics.trs import TrsVerifier
from .pipeline import stress_profile
from .thermostability import ThermostabilityVerifier

# the lightweight verifiers are stateless over read-only package data — build once and reuse,
# so a batch `rank` doesn't re-parse JSON per design. Geometry uses the sharpest reference on
# hand: the CSD metal–organic prior if it's been built, else the committed PDB reference.
REFERENCE = best_reference()
GEOMETRY = GeometryVerifier(REFERENCE)
BOND_VALENCE = BondValenceVerifier()
COORD_SYMMETRY = CoordinationSymmetryVerifier()  # nVECSUM: is the metal enclosed?
COORD_GEOMETRY = CoordinationGeometryVerifier()  # polyhedron shape vs ideal

# stages with a library verifier but no inline-available input (licence / prediction /
# scorer) — advertised so an agent knows the full stack and how to enable each.
NEEDS_INPUT = {
    "mogul": "a CSD licence (Mogul / CSD Python API)",
    "metalhawk": "MetalHawk geometry predictions (scripts/metalhawk_score.py — open, no licence; "
                 "EXPERIMENTAL: the learned ANN is confidently-OOD on de-novo designs — coord_geometry "
                 "is the analytic geometry oracle, see docs/experiments/2026-07-07-metalhawk-ood-designed-sites.md)",
    "cofold": "a co-fold prediction (scripts/chai_crosscheck, alphafold3_crosscheck, or allmetal3d_crosscheck)",
    "expression": "a sequence scorer (scripts/expression_score)",
    "thermostability": "an MD/Tm scorer (scripts/thermostability_score)",
    "motif_selectivity": "pass selectivity_metals (no GPU — MetalPDB donor-set enrichment: is this "
                         "donor set characteristic of the target metal in real proteins? an occupancy "
                         "prior, not a binding free energy)",
    "selectivity": "pass deep=True + selectivity_metals (MLIP metal-swap ΔE — does the target metal bind "
                   "best?). GATED + currently inert: every available MLIP backbone fails the Irving-Williams "
                   "series, so the tier defers rather than emit a meaningless metal ranking — see "
                   "docs/experiments/2026-07-13-mlip-cannot-rank-metals.md",
}


def as_dict(v: Verdict) -> dict:
    d = {"label": v.label, "score": round(v.score, 3), "trust": v.trust, "ood": v.ood, "reason": v.reason}
    if v.metrics:  # machine-readable numbers behind the verdict (σ, BVS, drift…) when the verifier exposes them
        d["metrics"] = v.metrics
    return d


MLIP_TIERS = ("mlip", "mlip_md", "trs")  # the deep tiers; `selectivity` joins them under selectivity_metals

# the full verifier stack in cost order — the unified `stack` view lists every tier with its
# status so the consensus is auditable, even tiers that didn't run on this input
STACK_ORDER = (
    "geometry", "bond_valence", "coord_symmetry", "coord_geometry", "precedent", "motif_selectivity",
    "metalhawk", "mogul", "mlip", "mlip_md", "trs", "selectivity", "cofold", "expression",
    "thermostability",
)


def stack(results: dict) -> list[dict]:
    """One entry per stack tier, in cost order: ran (with its verdict), skipped (backend
    absent), or needs_input (licence / scorer / co-fold) — the complete per-stage picture."""
    rows = []
    for stage in STACK_ORDER:
        if stage in results and "label" in results[stage]:
            rows.append({"stage": stage, "status": "ran", **results[stage]})
        elif stage in results:  # ran-but-skipped (e.g. MLIP with no backend)
            rows.append({"stage": stage, "status": "skipped", "detail": results[stage]["skipped"]})
        elif stage in NEEDS_INPUT:
            rows.append({"stage": stage, "status": "needs_input", "detail": NEEDS_INPUT[stage]})
        elif stage in MLIP_TIERS:  # only attempted with deep=True + a GPU backend
            rows.append({"stage": stage, "status": "needs_input", "detail": "pass deep=True (needs a GPU backend)"})
    return rows


AUTO = object()  # verify_structure default: build the MLIP backbone per call (a batch passes a shared one)

# tiers that exist only when their input is supplied: name → verifier factory
PROVIDED_TIERS: dict[str, Callable[[object], Verifier]] = {
    "cofold": CofoldCrossCheck,
    "metalhawk": MetalHawkVerifier,
    "expression": ExpressionVerifier,
    "mogul": MogulVerifier,
    "thermostability": ThermostabilityVerifier,
}

NO_MLIP_BACKEND = "no MLIP backend (install touchstone[mace])"


def default_verifiers(
    precedent: bool, precedent_search: Callable[..., object] | None, selectivity_metals: tuple[str, ...]
) -> dict[str, Verifier]:
    """The analytic tiers that run anywhere — no GPU, no licence."""
    tiers: dict[str, Verifier] = {
        "geometry": GEOMETRY,
        "bond_valence": BOND_VALENCE,
        "coord_symmetry": COORD_SYMMETRY,
        "coord_geometry": COORD_GEOMETRY,
    }
    if precedent:  # open MetalPDB coordination-motif precedent — default on (disable with precedent=False)
        tiers["precedent"] = PrecedentVerifier(precedent_search)
    if selectivity_metals:  # metal discrimination from observed occupancy — CPU-only, no deep needed
        tiers["motif_selectivity"] = MotifSelectivityVerifier(selectivity_metals)
    return tiers


def mlip_tiers(calc: object, selectivity_metals: tuple[str, ...]) -> dict[str, Verifier]:
    """The deep (GPU) tiers, sharing the one backbone (they protonate internally). Callers
    handle an absent backend — see `deep_skipped`."""
    tiers: dict[str, Verifier] = {
        "mlip": MLIPVerifier(calculator=calc),
        "mlip_md": MLIPDynamicsVerifier(calculator=calc),
        "trs": TrsVerifier(calculator=calc),  # preorganization: reorganization on unbinding
    }
    if selectivity_metals:  # MLIP metal-swap ΔE — does the target metal bind best?
        tiers["selectivity"] = MLIPSelectivityVerifier(calculator=calc, metals=selectivity_metals)
    return tiers


def deep_skipped(selectivity_metals: tuple[str, ...]) -> dict[str, dict]:
    """The deep tiers named as *skipped* when no backend is installed — skipped, not deferred,
    so an absent GPU never tanks the consensus."""
    names = MLIP_TIERS + (("selectivity",) if selectivity_metals else ())
    return {n: {"skipped": NO_MLIP_BACKEND} for n in names}


def run_tiers(verifiers: dict[str, Verifier], design: BinderDesign) -> tuple[dict[str, dict], list[str]]:
    """Run every tier → (results, the labels that count toward consensus). An unexpected
    per-verifier failure is recorded as skipped and excluded from the verdict, never raised."""
    results: dict[str, dict] = {}
    counted: list[str] = []
    for name, verifier in verifiers.items():
        try:
            verdict = verifier.verify(design)
        except Exception as e:  # unexpected per-verifier failure ⇒ skipped, not counted
            results[name] = {"skipped": f"{type(e).__name__}: {e}"}
            continue
        results[name] = as_dict(verdict)
        counted.append(verdict.label)
    return results, counted


def consensus_of(counted: list[str]) -> Literal["trust", "weak", "defer"]:
    """Aggregates the module docstring's rule over the tiers that actually ran."""
    if not counted or "defer" in counted:
        return "defer"
    return "trust" if all(label == "trust" for label in counted) else "weak"


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
    sequence: str = "",
    cofold_provider=None, metalhawk_scorer=None,
    precedent: bool = True, precedent_search=None,
    expression_scorer=None, mogul_analyse=None, thermostability_predictor=None,
    selectivity_metals=None, calc=AUTO,
) -> dict:
    """Verify a metal-coordination structure. Returns per-verifier verdicts, a `not_run`
    map of stages needing inputs, and a trust/weak/defer consensus. With `stress`, also
    re-verify the site under extreme-condition perturbations (acidic-leachate bond stretch,
    low-pH donor protonation) → a `stress` map {neutral/leachate/low_pH: verdict}: does it
    hold up in the real recovery process? `cofold_provider` (a design → predicted
    CoordinationSite callback over Chai-1 / AllMetal3D outputs) adds the independent co-fold
    cross-check tier; `metalhawk_scorer` (a design → MetalHawkPrediction over
    scripts/metalhawk_score.py output) adds the open geometry-distortion tier. The open MetalPDB
    **precedent** tier runs by default — disable with `precedent=False`, override the searcher with
    `precedent_search`. `expression_scorer`, `mogul_analyse` (licensed CSD), and
    `thermostability_predictor` each enable their opt-in tier — an absent backend never collapses the
    default consensus. The expression / thermostability scorers key by `sequence` (the site alone has
    none), so pass `sequence` to enable them.
    `selectivity_metals` (e.g. ("Ni2+","Cu2+","Co2+")) adds the MLIP metal-swap ΔE selectivity tier
    under `deep` — the physics discrimination geometry can't make (does the target bind best?). `calc`
    is an internal knob for batch callers (`rank_structures`) to share one MLIP backbone; leave it
    default."""
    site = coordination_site(structure, element_symbol(metal).upper(), metal, cutoff)
    design = BinderDesign(sequence, site, generator="external", generator_confidence=0.0, source=str(structure))

    metals = tuple(selectivity_metals or ())  # normalize once — both tier builders want a tuple
    verifiers = default_verifiers(precedent, precedent_search, metals)
    verifiers |= {name: PROVIDED_TIERS[name](p) for name, p in {
        "cofold": cofold_provider,  # independent predictor (Chai-1 / AllMetal3D) corroboration
        "metalhawk": metalhawk_scorer,  # independent ANN geometry-distortion oracle
        "expression": expression_scorer,  # sequence expressibility (ESM-2 pseudo-ppl + solubility)
        "mogul": mogul_analyse,  # licensed CSD Mogul geometry validation
        "thermostability": thermostability_predictor,  # whole-protein Tm (TemStaPro / DeepSTABp)
    }.items() if p is not None}

    skipped: dict[str, dict] = {}
    if deep:
        if calc is AUTO:  # single call ⇒ build per call; a batch hands in a shared backbone (or None)
            calc = mlip_backbone()
        if calc is None:  # no backend ⇒ skip, don't defer
            skipped = deep_skipped(metals)
        else:
            verifiers |= mlip_tiers(calc, metals)

    ran, counted = run_tiers(verifiers, design)
    results = skipped | ran
    consensus = consensus_of(counted)
    result = {
        "structure": str(structure),
        "metal": metal,
        "coordination_number": site.coordination_number,
        "donors": list(site.ligand_elems),
        "reference": REFERENCE.source,  # which geometry prior backed the z-score (CSD or PDB)
        "verifiers": results,
        "not_run": NEEDS_INPUT,
        "stack": stack(results),  # full per-tier breakdown (ran / skipped / needs_input), cost order
        "consensus": consensus,
    }
    if stress:  # robustness map: does the site hold its verdict across the operating envelope?
        # geometry-tier by design (independent of `deep`): the perturbations are geometric
        # (bond stretch, donor protonation), so the z-score is the natural judge — and running
        # the MLIP tier across every condition would multiply GPU cost for little extra signal.
        result["stress"] = {cond: as_dict(v) for cond, v in stress_profile(design, GEOMETRY).items()}
    return result
