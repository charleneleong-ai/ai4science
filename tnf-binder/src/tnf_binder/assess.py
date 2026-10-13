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

import re
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from pathlib import Path

from touchstone.core import Verdict
from touchstone.sampled import Sampled, stability

from tnf_binder.design import TnfDesign
from tnf_binder.structure import Complex
from tnf_binder.tables import design_stem, number
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


@dataclass(frozen=True)
class SampledAssessment:
    """One design assessed on every diffusion sample of its complex, not just the first.

    `assess` reads one structure, which is what the ranking used and what every structural
    count in this campaign's funnels was computed from. Re-read across all three samples, 10
    of 16 candidates change their clash-free verdict and interface area moves by a factor of
    two to five on the same design -- one reads 0 clashes over 658 A^2 on `model_0` and 9 over
    1203 A^2 on `model_2`.

    This reports that and changes nothing. A verdict recomputed on aggregated structures would
    reject a different set of designs, and swapping one selection for another after seeing
    both is not a correction.
    """

    name: str
    per_sample: tuple[Assessment, ...]

    @property
    def clashes(self) -> Sampled:
        # Named for the design, not the metric: `stability` reports which designs flipped.
        return Sampled.of(self.name, (a.metrics.clashes for a in self.per_sample if a.metrics))

    @property
    def interface_area(self) -> Sampled:
        return Sampled.of(self.name,
                          (a.metrics.area for a in self.per_sample if a.metrics))

    @property
    def clash_free_is_stable(self) -> bool:
        """False when the clash-free verdict depends on which sample was scored."""
        return self.clashes.is_stable(lambda c: c == 0)

    def caveat(self) -> str | None:
        """What to print beside a verdict that one more sample would have changed."""
        if self.clash_free_is_stable:
            return None
        return (f"{self.name}: clash-free verdict is sample-dependent "
                f"({[int(c) for c in self.clashes.values]}), "
                f"interface area {self.interface_area.spread:.0f} A^2 apart across samples")


def assess_samples(directory: Path, *, epitope: Collection[int] = POCKET,
                   divergent: Collection[int] | None = None) -> SampledAssessment:
    """Assess every `*_model_*.cif` under `directory` as one design."""
    paths = sorted(directory.glob("*_model_*.cif"))
    return SampledAssessment(directory.name,
                             tuple(assess(p, epitope=epitope, divergent=divergent)
                                   for p in paths))


def unstable_verdicts(assessments: Sequence[SampledAssessment]) -> str:
    """One line per design whose structural verdict rests on the sample that was scored.

    Run this before quoting a count from a structural gate. A funnel built on an unstable one
    reports the sample, not the designs.
    """
    report = stability([a.clashes for a in assessments], lambda c: c == 0)
    lines = [str(report)] + [c for a in assessments if (c := a.caveat())]
    return "\n".join(lines)


def geometry_across_samples(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """One geometry row per design, taken at its least favourable sample rather than its first.

    `score_designs` scores every structure; this used to keep only `_model_0`, which made each
    structural verdict a statement about one diffusion sample. Across the two slot-filling
    rounds, 10 of 16 designs change their clash-free verdict depending on which sample is read,
    and interface area moves by a factor of two to five on the same design.

    So clashes take the maximum, buried area the minimum, and `flag` the union of every flag
    any sample raised: a design that clashes on one structure clashes, and one whose interface
    slid off the epitope on one structure did that too. `sample_spread` carries the clash range
    so a caller can say which verdicts were close.

    **On the only data this has been tried against, it rejects everything.** All sixteen
    slot-filling candidates fail under it, at a clash ceiling of zero and of three, including
    both designs that were actually submitted: `helix0_R35H_E67H_E70H` reads epitope focus
    0.51, 0.18 and 0.30 against a 0.25 bar, and `helix0_E38H_E67H_E70H` reads 0.22, 0.22 and
    0.28. A gate that rejects everything discriminates exactly as much as one that accepts
    everything.

    That is not an argument for going back to reading `model_0`, which decided 10 of those 16
    verdicts by which structure was opened. It is the finding that neither aggregation works
    here: with three samples and this much variance, these structural gates cannot separate
    these designs, and a verdict either way is noise. What would change that is more samples
    before gating, or gating on the fraction of samples that pass rather than on any single
    one -- a different design from this, and unbuilt.
    """
    by_design: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        design = re.sub(r"_model_\d+$", "", design_stem(row["design"]))
        by_design.setdefault(design, []).append(row)
    merged = {}
    for design, samples in by_design.items():
        clashes = [int(number(r, "clashes")) for r in samples]
        flags = sorted({f for r in samples
                        for f in (r.get("flag") or "").strip().strip('"').split(",") if f})
        worst = max(samples, key=lambda r: int(number(r, "clashes")))
        focus = [number(r, "epitope_focus") for r in samples]
        merged[design] = {**worst,
                          "clashes": str(max(clashes)),
                          "interface_area": str(min(number(r, "interface_area") for r in samples)),
                          "epitope_focus": str(min(focus)),
                          "flag": ",".join(flags),
                          "n_samples": str(len(samples)),
                          "sample_spread": f"{min(clashes)}-{max(clashes)}",
                          "focus_spread": f"{min(focus):.2f}-{max(focus):.2f}"}
    return merged
