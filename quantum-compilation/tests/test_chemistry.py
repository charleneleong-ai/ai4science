"""The chemistry screen, against the numbers the real verifier reported for the same files.

Every test needs the 122 MB FCIDUMP, so the module skips without it. That keeps CI green on a
clone while still pinning the screen locally, which is where it is used.
"""

import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml
from quantum_compilation.chemistry import ChemistryScreen
from quantum_compilation.cost import FEMOCO_C0, MAX_HAMILTONIAN_ERROR_HA, lambda_eff

EXPERIMENTS_DIR = Path(__file__).parents[1] / "docs" / "experiments"
VERIFIERS = os.environ.get("VQCC_VERIFIERS", "")
VQCC_DATA = Path(os.environ.get("VQCC_DATA", Path.home() / "quantum-compilation-data"))
SUBMITTED = VQCC_DATA / "public/submission_v2/root/dfthcsos/submission.npz"

pytestmark = pytest.mark.skipif(
    not (VERIFIERS and Path(VERIFIERS).is_dir() and SUBMITTED.is_file()),
    reason="needs VQCC_VERIFIERS and the submitted factorization",
)


@pytest.fixture(scope="module")
def part() -> dict[str, object]:
    """What the official verifier recorded for submission 57834335's chemistry part."""
    record = yaml.safe_load((EXPERIMENTS_DIR / "57834335-parts.yaml").read_text())
    return record["parts"]["dfthcsos"]


@pytest.fixture(scope="module")
def screen() -> ChemistryScreen:
    return ChemistryScreen(VERIFIERS)


@pytest.fixture(scope="module")
def metrics(screen):
    return screen.evaluate(SUBMITTED)


@pytest.fixture(scope="module")
def resources(screen, part):
    """`T` and `Q` rebuilt from the shape and `Q_total` the parts record carries."""
    return screen.resources_from_record(part)


@pytest.fixture(scope="module")
def coarsened(screen, tmp_path_factory):
    """The submitted factorization with every angle re-rounded onto the `brot = 7` grid.

    A fresh fit at that precision would beat re-rounding, but not by the two orders of
    magnitude the result leaves to find, which is what makes this a wall and not an estimate.
    """
    with np.load(SUBMITTED, allow_pickle=False) as arrays:
        data = dict(arrays)
    brot = int(data["brot"]) - 1
    scaled = np.rint(np.asarray(data["angle_words"], dtype=np.float64) / 2.0)
    # Words are unsigned in [0, 2**brot); rounding the largest up lands one past the end.
    data["angle_words"] = np.clip(scaled, 0, 2**brot - 1).astype(data["angle_words"].dtype)
    data["brot"] = np.asarray(brot)
    path = tmp_path_factory.mktemp("coarse") / "submission.npz"
    np.savez(path, **data)
    return screen.evaluate(path)


class TestAgainstTheVerifier:
    """The metrics the screen reproduces from the words alone, read off the parts record."""

    @pytest.mark.parametrize(
        ("field", "recorded"),
        [("one_norm", "Lambda"), ("e_gap", "E_gap"), ("lambda_eff", "lambda_eff")],
    )
    def test_reproduces_the_recorded_metric(self, metrics, part, field, recorded):
        assert getattr(metrics, field) == pytest.approx(part[recorded], rel=1e-12)

    def test_reproduces_the_binding_residual(self, metrics, part):
        """The parts record keeps `max(error_1, error_2)`, which is the one the cap reads."""
        recorded = part["constraint_margins"]["hamiltonian_error_ha"]["value"]
        assert metrics.hamiltonian_error_ha == pytest.approx(recorded, rel=1e-12)

    def test_lambda_eff_follows_from_the_gap_and_the_one_norm(self, metrics):
        """The one relation worth asserting: the verifier's own `lambda_eff` formula."""
        assert metrics.lambda_eff == pytest.approx(
            lambda_eff(metrics.e_gap, metrics.one_norm), rel=1e-12
        )

    def test_the_cost_follows_once_the_counts_are_supplied(self, metrics, part):
        assert metrics.cost(part["T_total"], part["Q_total"]) == pytest.approx(
            part["our_measurement"]["cost"], rel=1e-9
        )


class TestDecomposition:
    """The split that picks the target: which term of `Lambda` is worth attacking."""

    def test_the_eigenvalue_masses_rebuild_the_one_norm(self, metrics):
        assert metrics.positive_mass + metrics.negative_mass + metrics.lambda_h2 == (
            pytest.approx(metrics.one_norm, rel=1e-12)
        )

    def test_the_gap_rebuilds_from_its_own_terms(self, metrics):
        assert metrics.c0 + 2 * metrics.negative_mass + metrics.half_sum_w2 == pytest.approx(
            metrics.e_gap, rel=1e-9
        )

    def test_the_two_body_one_norm_dominates(self, metrics):
        """Why the experiment targets `lambda_h2`: it is nearly all of `Lambda`."""
        assert metrics.lambda_h2 / metrics.one_norm > 0.9
        assert metrics.positive_mass < metrics.negative_mass < metrics.lambda_h2


class TestTheFloorOnTheScore:
    """`D`'s split into a movable and a fixed part -- the result that closed the refit."""

    def test_the_movable_mass_rebuilds_from_its_two_terms(self, metrics):
        """`D + c0 = 2P + 2*excess`, both non-negative -- the identity behind the floor.

        The two sides come from different quantities (`D` and `c0` against `P` and the
        `lambda_h2`/`W` excess), so agreement is a real check, and the excess being positive
        is what stops any factorization pushing `D` below `-c0`.
        """
        assert metrics.movable_mass == pytest.approx(
            2 * metrics.positive_mass + 2 * metrics.two_body_excess, rel=1e-14
        )
        assert metrics.two_body_excess > 0.0

    def test_c0_is_the_constant_cost_py_assumes(self, metrics):
        """The CI-green floor tests hardcode `FEMOCO_C0`; this is what pins it to the verifier."""
        assert metrics.c0 == pytest.approx(FEMOCO_C0, rel=1e-12)

    def test_the_whole_search_is_worth_three_percent(self, metrics):
        assert metrics.lambda_eff_headroom == pytest.approx(0.0312, abs=5e-4)


class TestRejection:
    """That a breach of each screenable constraint is named rather than passed over."""

    def test_the_submitted_factorization_is_accepted(self, metrics):
        assert metrics.rejected_for == ()
        assert {name for name in metrics.margins} == {
            "hamiltonian_error_ha",
            "e_gap_ha",
            "lcu_lambda",
        }

    @pytest.mark.parametrize(
        ("field", "bad", "named"),
        [
            ("error_2", 0.021, "hamiltonian_error_ha"),
            ("error_1", 0.021, "hamiltonian_error_ha"),
            ("e_gap", 3.9, "e_gap_ha"),
            ("one_norm", 48.0, "lcu_lambda"),
            ("lambda_eff", float("nan"), "gap_in_domain"),
        ],
    )
    def test_a_breached_constraint_is_named(self, metrics, field, bad, named):
        assert any(named in r for r in replace(metrics, **{field: bad}).rejected_for)


class TestShapeResources:
    """`T` and `Q` rebuilt from the shape, against what the verifier reported for that shape."""

    @pytest.mark.parametrize(("field", "recorded"), [("toffoli", "T_total"), ("qubits", "Q_total")])
    def test_reproduces_the_verifier(self, resources, part, field, recorded):
        assert getattr(resources, field) == part[recorded]

    def test_J_comes_out_of_the_recorded_total(self, resources, part):
        """The one term no report states on its own, so `Q` is only reproducible if it is right."""
        assert resources.angle_word_qubits == 5 * int(part["brot"])


class TestTheRotationPrecisionWall:
    """Why `brot` cannot be lowered: the submitted angles already sit on the 20 mHa cap."""

    def test_one_bit_coarser_breaches_the_cap_by_two_orders_of_magnitude(self, coarsened):
        assert coarsened.hamiltonian_error_ha > 100 * MAX_HAMILTONIAN_ERROR_HA
        assert any("hamiltonian_error_ha" in r for r in coarsened.rejected_for)
