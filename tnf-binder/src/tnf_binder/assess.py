"""Turning one predicted complex into a judged design — the step every consumer repeated.

Reading a complex and getting verdicts out of it is five steps that must happen in order and
with the same arguments: classify the chains, sum buried area per mature position, turn that
into interface metrics against a chosen epitope, build the design, then run the verifiers over
it. The CLI did this to render a table; the batch scorer did it to write a CSV row; both
spelled it out, which is how they came to disagree about which verifiers to run and which
epitope to measure against.

So the sequence lives here once and the callers only choose how to display it.

`epitope` is a parameter rather than a constant because the campaign measured focus against
both conserved patches and found no difference between them — a result that is only reportable
if the choice is explicit at the call site.
"""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path

from touchstone.core import Verdict

from tnf_binder.design import TnfDesign
from tnf_binder.structure import Complex
from tnf_binder.target import POCKET, TnfTarget
from tnf_binder.verifiers import CompositionVerifier, InterfaceVerifier, PkaSwitchVerifier
from tnf_binder.verifiers.interface import InterfaceMetrics


@dataclass(frozen=True)
class Assessment:
    """One complex, measured and judged. `design` is None when the chains could not be split."""

    name: str
    design: TnfDesign | None = None
    metrics: InterfaceMetrics | None = None
    verdicts: tuple[tuple[str, Verdict], ...] = ()

    @property
    def is_resolved(self) -> bool:
        """False when the binder could not be told apart from the target, so nothing was judged."""
        return self.design is not None

    @property
    def failures(self) -> tuple[str, ...]:
        """Why each verifier that JUDGED the design rejected it.

        Deferrals are excluded deliberately. `defer` means the verifier had no measurement to
        judge by, which is not a defect in the design, and the two must not be collapsed: a
        design with no pKa prediction is not a design whose pH switch was examined and found
        missing. Collapsing them put 'no pKa prediction available' in a column of reasons to
        reject, where it would have ranked designs down for a gap in our own inputs.
        """
        return tuple(v.reason for _, v in self.verdicts if v.label == "weak")

    @property
    def deferred(self) -> tuple[str, ...]:
        """What could not be judged at all, for the caller that wants to know it was asked."""
        return tuple(name for name, v in self.verdicts if v.label == "defer")


def assess(path: Path, *, epitope: Collection[int] = POCKET,
           divergent: Collection[int] | None = None) -> Assessment:
    """Measure and judge one predicted complex.

    `pka_switch` is included and will defer unless the caller has already attached histidine
    pKa measurements to the design — defer meaning "cannot judge", which is the honest verdict
    when no pKa was predicted, and distinct from judging the switch and finding it absent.
    """
    complex_ = Complex.from_file(path)
    if not complex_.is_resolved:
        return Assessment(path.stem)

    metrics = InterfaceMetrics.from_areas(
        complex_.buried_by_position(), epitope=epitope,
        divergent=TnfTarget().divergent if divergent is None else divergent,
        clashes=complex_.clashes())
    design = TnfDesign(path.stem, complex_.binder_sequence, structure=path,
                       metrics={"interface": metrics})
    verdicts = (("composition", CompositionVerifier().verify(design)),
                ("interface", InterfaceVerifier().verify(design)),
                ("pka_switch", PkaSwitchVerifier().verify(design)))
    return Assessment(path.stem, design, metrics, verdicts)
