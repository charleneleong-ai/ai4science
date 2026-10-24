from pathlib import Path

import pytest

tracking = pytest.importorskip("casmi.tracking", reason="needs the `tracking` extra")

VERDICTS = Path(__file__).parent.parent / "docs/experiments/casmi/verdicts.yaml"


def _sub(version: int, score: float) -> dict:
    return {"version": version, "score": score}


class TestAnnotate:
    """Best-so-far and gate G1 per submission."""

    def test_gate_needs_more_than_the_noise(self):
        out = tracking.annotate([_sub(1, 0.300), _sub(2, 0.304), _sub(3, 0.310), _sub(4, 0.200)], noise=0.005)
        assert [s["gate_g1"] for s in out] == [True, False, True, False]
        assert [s["best_before"] for s in out] == [None, 0.300, 0.304, 0.310]
        assert [s["best_so_far"] for s in out] == [0.300, 0.304, 0.310, 0.310]

    def test_recorded_verdicts_agree_with_the_gate(self):
        data = tracking.load(VERDICTS)
        for s in tracking.annotate(data["submissions"], data["noise"]):
            assert (s["verdict"] == "PROMOTE") == s["gate_g1"], f"{s['label']}: verdict disagrees with G1"


def test_narrative_keeps_sections_in_order(tmp_path):
    md = tmp_path / "report.md"
    md.write_text("## Highlights\nbest 0.349\n\n## Task\nrank structures\n")
    sections = tracking.narrative(md)
    assert list(sections) == ["Highlights", "Task"]
    assert sections["Highlights"] == "## Highlights\nbest 0.349\n"
