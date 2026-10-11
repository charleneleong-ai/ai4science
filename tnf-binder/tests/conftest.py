"""Fixtures shared across test modules.

MODEL_0 and MODEL_2 are verbatim Boltz-2 output from two diffusion samples of one BindCraft
trajectory (bc_traj0) folded against the TNF-alpha trimer: chain 0 is the binder, chains 1-3
the trimer. Keeping the real numbers means the artefacts they exhibit -- an asymmetric pair
matrix, and a global ipTM inflated by the target's own interfaces -- are tested as observed
rather than as invented.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tnf_binder.confidence import BoltzConfidence

MODEL_2 = {
    "iptm": 0.7146211266517639,
    "chains_ptm": {"0": 0.5384855270385742, "1": 0.9223341345787048,
                   "2": 0.9211609363555908, "3": 0.9266692996025085},
    "pair_chains_iptm": {
        "0": {"0": 0.5384855270385742, "1": 0.23248620331287384,
              "2": 0.21757104992866516, "3": 0.22476425766944885},
        "1": {"0": 0.20019255578517914, "1": 0.9223341345787048,
              "2": 0.8062747120857239, "3": 0.8299746513366699},
        "2": {"0": 0.17391160130500793, "1": 0.8206294775009155,
              "2": 0.9211609363555908, "3": 0.8071165084838867},
        "3": {"0": 0.14855808019638062, "1": 0.8250360488891602,
              "2": 0.8283896446228027, "3": 0.9266692996025085}},
}
MODEL_0 = {
    "iptm": 0.7707878947257996,
    "chains_ptm": {"0": 0.6596081852912903, "1": 0.9546123147010803,
                   "2": 0.9524770379066467, "3": 0.9551080465316772},
    "pair_chains_iptm": {
        "0": {"0": 0.6596081852912903, "1": 0.259974867105484,
              "2": 0.2480950802564621, "3": 0.24962498247623444},
        "1": {"0": 0.25198549032211304, "1": 0.9546123147010803,
              "2": 0.8763281106948853, "3": 0.8952808976173401},
        "2": {"0": 0.20570412278175354, "1": 0.8817487955093384,
              "2": 0.9524770379066467, "3": 0.8773713707923889},
        "3": {"0": 0.1804904341697693, "1": 0.8899710178375244,
              "2": 0.8938556909561157, "3": 0.9551080465316772}},
}


def write_sample(directory: Path, payload: dict[str, object], index: int = 0) -> Path:
    path = directory / f"confidence_model_{index}.json"
    path.write_text(json.dumps(payload))
    return path


@pytest.fixture
def predictions(tmp_path: Path) -> Path:
    """A predictions/ tree with one two-sample design, as Boltz writes it."""
    design = tmp_path / "bc_traj0"
    design.mkdir()
    for index, payload in ((0, MODEL_0), (2, MODEL_2)):
        (design / f"confidence_bc_traj0_model_{index}.json").write_text(json.dumps(payload))
    return tmp_path


@pytest.fixture
def sample(tmp_path: Path) -> BoltzConfidence:
    """One real diffusion sample: binder chain 0 against the TNF-alpha trimer."""
    return BoltzConfidence.from_file(write_sample(tmp_path, MODEL_2))
