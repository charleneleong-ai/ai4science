from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from touchstone.boltzgen import boltzgen_confidence

# three refold samples whose mean (0.55) sits 0.35 below the best (0.90) — the largest
# within-design spread measured over the 57-design campaign
SPREAD_SAMPLES = (0.375, 0.375, 0.90)


@pytest.fixture
def npz(tmp_path: Path) -> Callable[..., Path]:
    def write(name: str = "design", **arrays: object) -> Path:
        path = tmp_path / f"{name}.npz"
        np.savez(path, **arrays)
        return path

    return write


class TestAggregation:
    @pytest.mark.parametrize("aggregate, expected", [("mean", 0.55), ("best", 0.90)])
    def test_aggregate_picks_the_reduction_and_keeps_the_best_alongside(self, npz, aggregate, expected):
        conf = boltzgen_confidence(npz(design_to_target_iptm=SPREAD_SAMPLES), aggregate)
        assert conf["iptm"] == pytest.approx(expected)
        assert conf["iptm_best"] == pytest.approx(0.90)  # the gap stays visible either way

    def test_mean_is_the_default(self, npz):
        conf = boltzgen_confidence(npz(design_to_target_iptm=SPREAD_SAMPLES))
        assert conf["iptm"] == pytest.approx(0.55)

    def test_mean_and_best_rank_designs_differently(self, npz):
        # the measured reorder: one lucky sample puts `spiky` first on best-of-3 and last on the mean
        spiky = npz("spiky", design_to_target_iptm=(0.90, 0.10, 0.10))
        steady = npz("steady", design_to_target_iptm=(0.50, 0.50, 0.50))

        def order(aggregate: str) -> list[str]:
            scored = [(p.stem, boltzgen_confidence(p, aggregate)["iptm"]) for p in (spiky, steady)]
            return [name for name, _ in sorted(scored, key=lambda s: s[1], reverse=True)]

        assert order("mean") == ["steady", "spiky"]
        assert order("best") == ["spiky", "steady"]

    def test_scalar_confidence_is_treated_as_one_sample(self, npz):
        assert boltzgen_confidence(npz(iptm=np.float32(0.7)))["iptm"] == pytest.approx(0.7, abs=1e-6)


class TestKeyFallback:
    @pytest.mark.parametrize(
        "signal, arrays, expected",
        [
            ("iptm", {"design_to_target_iptm": (0.4,), "iptm": (0.9,)}, 0.4),  # preferred key wins
            ("iptm", {"ligand_iptm": (0.4,), "iptm": (0.9,)}, 0.4),
            ("iptm", {"iptm": (0.9,)}, 0.9),  # last resort
            ("iptm", {"complex_plddt": (0.8,)}, None),  # no iPTM under any name
            ("ptm", {"design_ptm": (0.3,), "ptm": (0.9,)}, 0.3),
            ("ptm", {"ptm": (0.9,)}, 0.9),
        ],
        ids=["iptm_primary", "iptm_ligand", "iptm_last_resort", "iptm_absent", "ptm_primary", "ptm_fallback"],
    )
    def test_falls_back_across_key_names(self, npz, signal, arrays, expected):
        assert boltzgen_confidence(npz(**arrays))[signal] == expected

    def test_zero_confidence_is_a_value_not_a_missing_key(self, npz):
        # a genuine 0.0 must not fall through to the next key name — that would report a
        # failed fold as a confident one
        conf = boltzgen_confidence(npz(design_to_target_iptm=(0.0, 0.0), ligand_iptm=(0.8,)))
        assert conf["iptm"] == 0.0


def test_missing_file_is_none(tmp_path):
    assert boltzgen_confidence(tmp_path / "nope.npz") is None
