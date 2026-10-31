import kickmix as km
import numpy as np
import pytest
from quantum_compilation.submission import NPZ_ARRAYS, audit

TINY = "CCX q0 q1 q2"


def write_circuit(path, suffix=".kmx"):
    path.parent.mkdir(parents=True, exist_ok=True)
    km.Circuit(TINY).to_file(str(path.with_suffix(suffix)))


def write_npz(path, arrays=NPZ_ARRAYS):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **{name: np.zeros(2) for name in arrays})


@pytest.fixture
def valid(tmp_path):
    """A structurally complete submission: two task circuits, three angles, a factorization."""
    write_circuit(tmp_path / "invsqrt")
    write_circuit(tmp_path / "dqi")
    for index in range(3):
        write_circuit(tmp_path / "dfthcsos" / f"angle_k_{index}")
    write_npz(tmp_path / "dfthcsos" / "submission.npz")
    return tmp_path


class TestValidSubmission:
    """A complete directory audits clean and reports what it found."""

    def test_reports_no_problems(self, valid):
        assert audit(valid).problems == ()
        assert audit(valid).ok

    def test_finds_every_circuit_in_index_order(self, valid):
        report = audit(valid)
        assert [p.stem for p in report.angle_circuits] == [
            "angle_k_0",
            "angle_k_1",
            "angle_k_2",
        ]
        assert set(report.circuits) == {
            "invsqrt",
            "dqi",
            "angle_k_0",
            "angle_k_1",
            "angle_k_2",
        }

    def test_accepts_the_binary_suffix_too(self, tmp_path, valid):
        (valid / "invsqrt.kmx").rename(valid / "invsqrt.kmb")
        assert audit(valid).invsqrt.name == "invsqrt.kmb"

    def test_orders_angle_circuits_numerically_not_lexically(self, valid):
        # angle_k_10 sorts before angle_k_2 as text; the contract is "numbered from 0".
        for index in (3, 10):
            write_circuit(valid / "dfthcsos" / f"angle_k_{index}")
        report = audit(valid)
        assert [p.stem for p in report.angle_circuits][-2:] == [
            "angle_k_3",
            "angle_k_10",
        ]
        assert any("missing index(es)" in p for p in report.problems)


class TestLayoutProblems:
    @pytest.mark.parametrize("stem", ["invsqrt", "dqi"])
    def test_a_missing_task_circuit_is_reported(self, valid, stem):
        (valid / f"{stem}.kmx").unlink()
        assert any(stem in problem for problem in audit(valid).problems)

    def test_a_gap_in_the_angle_numbering_is_reported(self, valid):
        (valid / "dfthcsos" / "angle_k_1.kmx").unlink()
        assert any("missing index(es) [1]" in p for p in audit(valid).problems)

    def test_angles_not_starting_at_zero_are_reported(self, valid):
        (valid / "dfthcsos" / "angle_k_0.kmx").unlink()
        assert any("missing index(es) [0]" in p for p in audit(valid).problems)

    def test_an_empty_chemistry_directory_is_reported(self, tmp_path):
        write_circuit(tmp_path / "invsqrt")
        write_circuit(tmp_path / "dqi")
        (tmp_path / "dfthcsos").mkdir()
        assert any("no angle_k_<n> circuits" in p for p in audit(tmp_path).problems)

    def test_a_missing_factorization_is_reported(self, valid):
        (valid / "dfthcsos" / "submission.npz").unlink()
        assert any("submission.npz" in p for p in audit(valid).problems)

    def test_a_factorization_missing_an_array_is_reported(self, valid):
        write_npz(valid / "dfthcsos" / "submission.npz", arrays=NPZ_ARRAYS[:-1])
        problems = audit(valid).problems
        assert any("tree_edges" in p for p in problems)

    def test_an_unreadable_factorization_is_reported_not_raised(self, valid):
        (valid / "dfthcsos" / "submission.npz").write_text("not an npz")
        assert any("could not be read as an npz" in p for p in audit(valid).problems)

    def test_a_missing_directory_is_a_single_problem(self, tmp_path):
        report = audit(tmp_path / "absent")
        assert len(report.problems) == 1
        assert "not a directory" in report.problems[0]

    def test_an_ambiguous_task_circuit_counts_as_missing(self, valid):
        # Both suffixes present: the metric takes one file per task, so this is not a
        # submission the leaderboard would accept and the audit must not pick arbitrarily.
        write_circuit(valid / "invsqrt", suffix=".kmb")
        assert audit(valid).invsqrt is None
