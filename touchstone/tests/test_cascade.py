import numpy as np

from touchstone import BinderDesign, cascade
from touchstone.core import CoordinationSite, Verdict


def _design(seq: str) -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
    return BinderDesign(seq, site, generator="test", generator_confidence=0.5)


class _Gate:
    """Cheap tier: off-manifold (defer) for any sequence not in `passes`."""

    def __init__(self, passes: set[str]):
        self.passes = passes

    def verify(self, design: BinderDesign) -> Verdict:
        ok = design.sequence in self.passes
        return Verdict(1.0 if ok else 0.0, trust=ok, ood=not ok, reason="gate")


class _WeakGate:
    """Cheap tier that returns weak (not trusted, but not off-manifold)."""

    def verify(self, design: BinderDesign) -> Verdict:
        return Verdict(0.3, trust=False, ood=False, reason="weak")


class _Spy:
    """Expensive tier: records which designs actually reached it."""

    def __init__(self):
        self.seen: list[str] = []

    def verify(self, design: BinderDesign) -> Verdict:
        self.seen.append(design.sequence)
        return Verdict(1.0, trust=True, ood=False, reason="spy")


class TestCascade:
    def test_expensive_tier_only_runs_on_survivors(self):
        designs = [_design(s) for s in ("A", "B", "C")]
        spy = _Spy()
        results = cascade(designs, [("gate", _Gate({"B"})), ("expensive", spy)])
        assert spy.seen == ["B"]  # the rejected designs never reached the expensive tier
        by_seq = {r.design.sequence: r for r in results}
        assert by_seq["A"].survived is False and by_seq["A"].dropped_at == "gate"
        assert len(by_seq["A"].verdicts) == 1  # short-circuited at the gate
        assert by_seq["B"].survived is True and by_seq["B"].dropped_at is None
        assert len(by_seq["B"].verdicts) == 2

    def test_advances_predicate_controls_strictness(self):
        designs = [_design("X")]
        lenient = _Spy()
        cascade(designs, [("gate", _WeakGate()), ("expensive", lenient)])
        assert lenient.seen == ["X"]  # weak advances under the default (drop only off-manifold)

        strict = _Spy()
        cascade(designs, [("gate", _WeakGate()), ("expensive", strict)], advances=lambda v: v.trust)
        assert strict.seen == []  # weak dropped under aggressive (trust-only) gating
