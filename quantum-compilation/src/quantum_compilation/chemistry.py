"""Screening the chemistry task's constraints without building a single circuit.

A full `dfthcsos` verification takes about a minute because it constructs the ANGLE circuits
and simulates them to recover the angle words. But every constraint that rejects a candidate
-- both Frobenius residuals, the one-norm floor, the energy gap -- and the first factor of the
score, `lambda_eff`, are functions of the *words* alone. Feeding the submitted words straight
into the verifier's own evaluator gives all of them in well under a second, which is what makes
a search over factorizations affordable at all.

Two honest limits, both load-bearing:

* **This is a screen, not a verification.** The real verifier scores the words it *recovers by
  simulating the circuits*, not the words the npz carries. On the submitted plateau
  factorization the two agree exactly -- to relative error 0.0 -- on `error_2`, `Lambda`,
  `E_gap` and `lambda_eff`, which is every metric that submission's record preserves. Its
  `error_1` was not recorded, so the drift between screened and recovered words is *unmeasured*
  on this instance; `tests/test_chemistry.py` pins it on the sample submission instead, whose
  committed verifier output carries both residuals separately.
* **It cannot produce a cost.** The score is `lambda_eff * T * Q`, and `T` counts Toffolis as a
  simulated per-shot mean. No screen can know it, so this one reports `lambda_eff` and the
  margins and stops there, the way `submission.Resources` reports counts and calls them a bound.

Unlike `verify`, this imports the official verifier in-process rather than running it under
`python -I`. The isolation that buys is about the *submission*, which is data: here the npz is
read with `allow_pickle=False` and nothing else from the candidate is executed. The verifier
package itself is the competition's own scorer, trusted to the same degree whether it is
imported or spawned.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .cost import DFTHCSOS_CONSTRAINTS, Constraint, dfthcsos_cost
from .submission import CHEMISTRY_DIR, NPZ_SCREEN_ARRAYS, SubmissionAudit
from .verify import verifier_script

# The cap that covers both residuals, so which metrics it reads stays stated in one place.
HAMILTONIAN_ERROR = next(c for c in DFTHCSOS_CONSTRAINTS if c.name == "hamiltonian_error_ha")


@dataclass(frozen=True)
class ChemistryMetrics:
    """What the verifier reports for a factorization, and the terms behind it.

    Only independent quantities are stored. `Lambda`'s split into eigenvalue mass, and the
    gap's into its own terms, are algebra on these -- see `negative_mass` -- so there is no
    second copy of a number to fall out of step with the first.
    """

    error_1: float
    error_2: float
    one_norm: float
    e_gap: float
    lambda_eff: float
    lambda_h2: float
    half_sum_w2: float
    c0: float

    @property
    def gap_in_domain(self) -> bool:
        """The verifier returns `lambda_eff` as NaN exactly when the gap leaves its domain."""
        return not math.isnan(self.lambda_eff)

    @property
    def negative_mass(self) -> float:
        """Negative eigenvalue mass of `h1_prime`, inverted out of `E_gap`.

        `E_gap = c0 + 2*negative_mass + half_sum_w2` is the verifier's own `E_SOS`, rearranged,
        so this recovers the term without repeating its eigen-decomposition.
        """
        return (self.e_gap - self.c0 - self.half_sum_w2) / 2.0

    @property
    def positive_mass(self) -> float:
        """Positive eigenvalue mass, from `Lambda = positive + negative + lambda_h2`."""
        return self.one_norm - self.lambda_h2 - self.negative_mass

    @property
    def gap_product(self) -> float:
        """`D = 2*Lambda - E_gap`. With `E_gap` pinned at its floor, this is what the score
        turns on, and `lambda_eff = sqrt(E_gap * D)` is a downward parabola in `E_gap`."""
        return 2.0 * self.one_norm - self.e_gap

    @property
    def report(self) -> dict[str, dict[str, float]]:
        """These metrics in the shape the chemistry verifier reports them.

        So `cost`'s constraint table reads a screened candidate through exactly the same
        `Constraint.paths` it reads a verified one through, instead of this module keeping a
        second copy of which metric each limit bounds.
        """
        return {
            "metrics": {
                "error_1": self.error_1,
                "error_2": self.error_2,
                "E_gap": self.e_gap,
                "Lambda": self.one_norm,
            }
        }

    @property
    def screenable(self) -> tuple[Constraint, ...]:
        """The constraints this screen can speak to at all.

        `ANGLE_WIDTH` drops out structurally rather than by name: it reads a `resources`
        path, and a screen that never builds the circuits has no resources to report.
        """
        metrics = self.report["metrics"]
        return tuple(
            c
            for c in DFTHCSOS_CONSTRAINTS
            if all(section == "metrics" and key in metrics for section, key in c.paths)
        )

    @property
    def hamiltonian_error_ha(self) -> float:
        """The binding residual: both are capped at the same 20 mHa, so the larger one binds."""
        return HAMILTONIAN_ERROR.value(self.report)

    @property
    def margins(self) -> dict[str, dict[str, Any]]:
        """Headroom against every constraint a screen can speak to."""
        report = self.report
        return {c.name: c.margin(report) for c in self.screenable}

    @property
    def rejected_for(self) -> tuple[str, ...]:
        """Every screenable constraint this candidate fails, by `verify`'s own table."""
        report = self.report
        reasons = [
            f"{c.name} {c.value(report):.6f} against {c.limit:g} ({c.sense})"
            for c in self.screenable
            if c.violated(c.value(report))
        ]
        if not self.gap_in_domain:
            reasons.append("gap_in_domain is false")
        return tuple(reasons)

    def cost(self, toffolis: float, qubits: int) -> float:
        """The score, once a verifier run has supplied the two counts a screen cannot know."""
        return dfthcsos_cost(qubits, toffolis, self.e_gap, self.one_norm)


class ChemistryScreen:
    """The FeMoCo reference, loaded once, against which many factorizations are screened.

    Loading the 122 MB FCIDUMP costs about four seconds and is the whole reason this is a
    class: a search evaluates thousands of candidates against one reference.
    """

    def __init__(self, verifiers_root: Path | str) -> None:
        # Goes through `verify`'s table so the path error is the one the rest of the CLI
        # gives, and so an empty `dfthcsos/` directory fails here rather than as an
        # ImportError further in.
        self.package_root = verifier_script(verifiers_root, CHEMISTRY_DIR).parent.parent
        if str(self.package_root) not in sys.path:
            sys.path.insert(0, str(self.package_root))

    def evaluate(self, factorization: Path | str) -> ChemistryMetrics:
        # Deferred: these resolve only once `__init__` has put the verifier on `sys.path`.
        from dfthcsos.scoreinfra.dfthc_math import evaluate_dfthc_angle_words, lambda_h2
        from dfthcsos.verify import load_inputs

        path = Path(factorization)
        if problems := SubmissionAudit.npz_problems(path, NPZ_SCREEN_ARRAYS):
            raise ValueError("; ".join(problems))

        inputs = load_inputs(path)
        with np.load(path, allow_pickle=False) as arrays:
            words, edges = arrays["angle_words"], arrays["tree_edges"]
            brot = int(arrays["brot"])
        reported = evaluate_dfthc_angle_words(
            words,
            inputs.w,
            inputs.W,
            inputs.S,
            inputs.beta_1,
            inputs.ref,
            brot=brot,
            tree_edges=edges,
            energy_upper=inputs.energy_upper,
        )
        outer_w = np.asarray(inputs.W, dtype=np.float64)
        ref = inputs.ref
        return ChemistryMetrics(
            error_1=float(reported["error_1"]),
            error_2=float(reported["error_2"]),
            one_norm=float(reported["Lambda"]),
            e_gap=float(reported["E_gap"]),
            lambda_eff=float(reported["lambda_eff"]),
            lambda_h2=float(lambda_h2(inputs.w, outer_w)),
            half_sum_w2=float(0.5 * np.sum(outer_w**2)),
            c0=float(
                inputs.energy_upper - ref["ecore"] + ref["electrons"] * float(inputs.beta_1)
            ),
        )
