"""Tests for touchstone.boltz — reading Boltz confidence with the binder/target split.

MODEL_0 and MODEL_2 are verbatim Boltz output from two diffusion samples of one design folded
against a homotrimeric target: chain 0 is the binder, chains 1-3 the trimer. Real numbers,
because the artefacts they exhibit are the point -- an asymmetric pair matrix, and a global
ipTM inflated by the target folding against itself.
"""

from __future__ import annotations

import json

import pytest

from touchstone.boltz import BoltzConfidence, DesignConfidence, rank_designs

MODEL_2 = {
    "iptm": 0.7146211266517639,
    "chains_ptm": {
        "0": 0.5384855270385742,
        "1": 0.9223341345787048,
        "2": 0.9211609363555908,
        "3": 0.9266692996025085,
    },
    "pair_chains_iptm": {
        "0": {
            "0": 0.5384855270385742,
            "1": 0.23248620331287384,
            "2": 0.21757104992866516,
            "3": 0.22476425766944885,
        },
        "1": {
            "0": 0.20019255578517914,
            "1": 0.9223341345787048,
            "2": 0.8062747120857239,
            "3": 0.8299746513366699,
        },
        "2": {
            "0": 0.17391160130500793,
            "1": 0.8206294775009155,
            "2": 0.9211609363555908,
            "3": 0.8071165084838867,
        },
        "3": {
            "0": 0.14855808019638062,
            "1": 0.8250360488891602,
            "2": 0.8283896446228027,
            "3": 0.9266692996025085,
        },
    },
}
MODEL_0 = {
    "iptm": 0.7707878947257996,
    "chains_ptm": {
        "0": 0.6596081852912903,
        "1": 0.9546123147010803,
        "2": 0.9524770379066467,
        "3": 0.9551080465316772,
    },
    "pair_chains_iptm": {
        "0": {
            "0": 0.6596081852912903,
            "1": 0.259974867105484,
            "2": 0.2480950802564621,
            "3": 0.24962498247623444,
        },
        "1": {
            "0": 0.25198549032211304,
            "1": 0.9546123147010803,
            "2": 0.8763281106948853,
            "3": 0.8952808976173401,
        },
        "2": {
            "0": 0.20570412278175354,
            "1": 0.8817487955093384,
            "2": 0.9524770379066467,
            "3": 0.8773713707923889,
        },
        "3": {
            "0": 0.1804904341697693,
            "1": 0.8899710178375244,
            "2": 0.8938556909561157,
            "3": 0.9551080465316772,
        },
    },
}


def write_sample(directory, payload, index=0):
    path = directory / f"confidence_model_{index}.json"
    path.write_text(json.dumps(payload))
    return path


@pytest.fixture
def sample(tmp_path) -> BoltzConfidence:
    """One real diffusion sample: binder chain 0 against the trimer."""
    return BoltzConfidence.from_file(write_sample(tmp_path, MODEL_2))


@pytest.fixture
def predictions(tmp_path):
    """A predictions/ tree with one two-sample design, as Boltz writes it."""
    design = tmp_path / "design_a"
    design.mkdir()
    for index, payload in ((0, MODEL_0), (2, MODEL_2)):
        (design / f"confidence_design_a_model_{index}.json").write_text(
            json.dumps(payload)
        )
    return tmp_path


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
        payload = {
            "iptm": 0.4,
            "chains_ptm": {"0": 0.9, "1": 0.9},
            "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.4}, "1": {"0": 0.4, "1": 0.9}},
        }
        sample = BoltzConfidence.from_file(write_sample(tmp_path, payload))
        assert sample.target_internal == 0.0
        assert sample.binder_to_target == pytest.approx(0.4)


class TestDesignConfidence:
    """Aggregation across diffusion samples."""

    def test_ranks_on_mean_not_best_of_n(self, predictions):
        design = DesignConfidence.from_dir(predictions / "design_a")
        assert len(design.samples) == 2
        assert design.binder_to_target == pytest.approx(0.216, abs=1e-3)
        assert design.binder_to_target_best == pytest.approx(0.233, abs=1e-3)
        assert design.binder_to_target < design.binder_to_target_best

    def test_inflation_is_the_gap_global_iptm_hides(self, predictions):
        design = DesignConfidence.from_dir(predictions / "design_a")
        assert design.global_iptm == pytest.approx(0.743, abs=1e-3)
        assert design.inflation == pytest.approx(
            design.global_iptm - design.binder_to_target
        )
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
        payload = {
            **MODEL_2,
            "pair_chains_iptm": {
                **MODEL_2["pair_chains_iptm"],
                "0": {"0": 0.5, "1": 0.05, "2": 0.05, "3": 0.05},
            },
        }
        write_sample(weak, payload)
        (predictions / "no_confidence_yet").mkdir()

        ranked = rank_designs(predictions)
        assert [d.name for d in ranked] == ["design_a", "weak_design"]
        assert ranked[1].binder_to_target < 0.15
