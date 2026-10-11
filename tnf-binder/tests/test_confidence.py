"""Boltz-2 confidence reading, anchored to real output from the campaign.

The MODEL_0 / MODEL_2 payloads and the `predictions` tree live in conftest.py, which documents
where they came from.
"""

from __future__ import annotations

import pytest

from tnf_binder.confidence import BoltzConfidence, DesignConfidence, rank_designs
from conftest import MODEL_2, write_sample


class TestBoltzConfidence:
    """One diffusion sample, with the binder/target split applied."""

    def test_global_iptm_overstates_the_designed_interface(self, sample):
        assert sample.global_iptm == pytest.approx(0.715, abs=1e-3)
        assert sample.binder_to_target == pytest.approx(0.200, abs=1e-3)
        assert sample.target_internal == pytest.approx(0.820, abs=1e-3)

    def test_pair_is_symmetrised(self, sample):
        """Boltz's matrix is asymmetric; neither direction is privileged."""
        assert sample.pair_iptm["0"]["3"] != sample.pair_iptm["3"]["0"]
        assert sample.pair("0", "3") == sample.pair("3", "0")
        assert sample.pair("0", "3") == pytest.approx(0.187, abs=1e-3)

    def test_target_chains_exclude_the_binder(self, sample):
        assert sample.target_chains == ("1", "2", "3")
        assert sample.binder_ptm == pytest.approx(0.538, abs=1e-3)

    def test_monomer_target_has_no_internal_interface(self, tmp_path):
        """A single-chain target cannot inflate the global score; nothing should divide by zero."""
        payload = {"iptm": 0.4, "chains_ptm": {"0": 0.9, "1": 0.9},
                   "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.4}, "1": {"0": 0.4, "1": 0.9}}}
        sample = BoltzConfidence.from_file(write_sample(tmp_path, payload))
        assert sample.target_internal == 0.0
        assert sample.binder_to_target == pytest.approx(0.4)


class TestDesignConfidence:
    """Aggregation across diffusion samples."""

    def test_ranks_on_mean_not_best_of_n(self, predictions):
        design = DesignConfidence.from_dir(predictions / "bc_traj0")
        assert len(design.samples) == 2
        assert design.binder_to_target == pytest.approx(0.216, abs=1e-3)
        assert design.binder_to_target_best == pytest.approx(0.233, abs=1e-3)
        assert design.binder_to_target < design.binder_to_target_best

    def test_inflation_is_the_gap_global_iptm_hides(self, predictions):
        design = DesignConfidence.from_dir(predictions / "bc_traj0")
        assert design.global_iptm == pytest.approx(0.743, abs=1e-3)
        assert design.inflation == pytest.approx(design.global_iptm - design.binder_to_target)
        assert design.inflation > 0.5

    def test_empty_directory_is_falsy_and_scores_zero(self, tmp_path):
        design = DesignConfidence.from_dir(tmp_path)
        assert not design
        assert design.binder_to_target == 0.0
        assert design.inflation == 0.0


class TestRankDesigns:
    """Scanning a predictions tree."""

    def test_orders_by_mean_interface_and_skips_empty_dirs(self, predictions):
        weak = predictions / "weak_design"
        weak.mkdir()
        payload = {**MODEL_2, "pair_chains_iptm": {
            **MODEL_2["pair_chains_iptm"],
            "0": {"0": 0.5, "1": 0.05, "2": 0.05, "3": 0.05}}}
        write_sample(weak, payload)
        (predictions / "no_confidence_yet").mkdir()

        ranked = rank_designs(predictions)
        assert [d.name for d in ranked] == ["bc_traj0", "weak_design"]
        assert ranked[1].binder_to_target < 0.15
