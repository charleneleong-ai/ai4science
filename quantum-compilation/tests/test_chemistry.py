"""The chemistry screen, against the numbers the real verifier reported for the same files.

Every test needs the 122 MB FCIDUMP, so the module skips without it. That keeps CI green on a
clone while still pinning the screen locally, which is where it is used.
"""

import os
from dataclasses import replace
from pathlib import Path

import pytest
import yaml
from quantum_compilation.chemistry import ChemistryScreen
from quantum_compilation.cost import lambda_eff

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
def metrics(part):
    return ChemistryScreen(VERIFIERS).evaluate(SUBMITTED)


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
