"""A submission's layout, and the measured content of the circuits in it."""

import kickmix as km
import numpy as np
import pytest
from quantum_compilation.submission import NPZ_ARRAYS, Resources, SubmissionAudit

TINY = "CCX q0 q1 q2"

# Two Toffoli-likes (a CCX and a CSWAP), four T-level rotations, an S gate that is Clifford,
# and one level-4 rotation that Clifford+T cannot express exactly.
MIXED = """
    X q0
    CCX q0 q1 q2
    CSWAP q0 q3 q4
    Z_POW q0 0.25
    Z_POW q1 -0.25
    Z_POW q2 0.75
    Z_POW q3 0.25
    Z_POW q4 0.125
    Z_POW q5 0.5
"""


def write_circuit(path, suffix=".kmx"):
    path.parent.mkdir(parents=True, exist_ok=True)
    km.Circuit(TINY).to_file(str(path.with_suffix(suffix)))


def write_npz(path, arrays=NPZ_ARRAYS):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **{name: np.zeros(2) for name in arrays})


def complains(root, needle: str) -> bool:
    """Whether auditing `root` reports a problem mentioning `needle`."""
    return any(needle in problem for problem in SubmissionAudit.from_root(root).problems)


@pytest.fixture
def valid(tmp_path):
    """A structurally complete submission: two task circuits, three angles, a factorization."""
    write_circuit(tmp_path / "invsqrt")
    write_circuit(tmp_path / "dqi")
    for index in range(3):
        write_circuit(tmp_path / "dfthcsos" / f"angle_k_{index}")
    write_npz(tmp_path / "dfthcsos" / "submission.npz")
    return tmp_path


@pytest.fixture(scope="module")
def mixed_circuit() -> km.Circuit:
    return km.Circuit(MIXED)


@pytest.fixture(scope="module")
def mixed(mixed_circuit) -> Resources:
    return Resources.from_circuit(mixed_circuit)


class TestValidSubmission:
    """A complete directory audits clean and reports what it found."""

    def test_reports_no_problems(self, valid):
        report = SubmissionAudit.from_root(valid)
        assert report.problems == ()
        assert report.ok

    def test_finds_every_circuit_in_index_order(self, valid):
        report = SubmissionAudit.from_root(valid)
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

    def test_accepts_the_binary_suffix_too(self, valid):
        (valid / "invsqrt.kmx").rename(valid / "invsqrt.kmb")
        assert SubmissionAudit.from_root(valid).invsqrt.name == "invsqrt.kmb"

    def test_orders_angle_circuits_numerically_not_lexically(self, valid):
        # angle_k_10 sorts before angle_k_2 as text; the contract is "numbered from 0".
        for index in (3, 10):
            write_circuit(valid / "dfthcsos" / f"angle_k_{index}")
        report = SubmissionAudit.from_root(valid)
        assert [p.stem for p in report.angle_circuits][-2:] == ["angle_k_3", "angle_k_10"]
        assert any("missing index(es)" in p for p in report.problems)


class TestLayoutProblems:
    @pytest.mark.parametrize("stem", ["invsqrt", "dqi"])
    def test_a_missing_task_circuit_is_reported(self, valid, stem):
        (valid / f"{stem}.kmx").unlink()
        assert complains(valid, stem)

    @pytest.mark.parametrize("removed", [0, 1])
    def test_a_gap_in_the_angle_numbering_is_reported(self, valid, removed):
        # Numbering must run from 0 with no gaps, so a missing first circuit is the same
        # failure as a missing middle one.
        (valid / "dfthcsos" / f"angle_k_{removed}.kmx").unlink()
        assert complains(valid, f"missing index(es) [{removed}]")

    def test_an_empty_chemistry_directory_is_reported(self, tmp_path):
        write_circuit(tmp_path / "invsqrt")
        write_circuit(tmp_path / "dqi")
        (tmp_path / "dfthcsos").mkdir()
        assert complains(tmp_path, "no angle_k_<n> circuits")

    def test_a_missing_factorization_is_reported(self, valid):
        (valid / "dfthcsos" / "submission.npz").unlink()
        assert complains(valid, "submission.npz")

    def test_a_factorization_missing_an_array_is_reported(self, valid):
        write_npz(valid / "dfthcsos" / "submission.npz", arrays=NPZ_ARRAYS[:-1])
        assert complains(valid, "tree_edges")

    def test_an_unreadable_factorization_is_reported_not_raised(self, valid):
        (valid / "dfthcsos" / "submission.npz").write_text("not an npz")
        assert complains(valid, "could not be read as an npz")

    def test_a_missing_directory_is_a_single_problem(self, tmp_path):
        report = SubmissionAudit.from_root(tmp_path / "absent")
        assert len(report.problems) == 1
        assert "not a directory" in report.problems[0]

    def test_an_ambiguous_task_circuit_counts_as_missing(self, valid):
        # Both suffixes present: the metric takes one file per task, so this is not a
        # submission the leaderboard would accept and the audit must not pick arbitrarily.
        write_circuit(valid / "invsqrt", suffix=".kmb")
        assert SubmissionAudit.from_root(valid).invsqrt is None


class TestResourceExtraction:
    """Width and non-Clifford content, read off the circuit by kickmix.

    `num_toffoli_like` is the static count -- what the verifiers report as
    `max_toffoli_count`. The scored figure is the simulated per-shot mean, which
    measurement-based uncomputation puts below it, so these counts bound a candidate's cost
    from above and never give it.
    """

    def test_counts_toffoli_likes_including_cswap(self, mixed):
        assert mixed.num_toffoli_like == 2

    def test_separates_t_level_from_finer_rotations(self, mixed):
        # The S gate is Clifford, so it is neither a T nor a fine rotation.
        assert mixed.num_t_like == 4
        assert mixed.num_fine_rotations == 1

    def test_counts_fine_rotations_across_clifford_levels(self):
        # Levels 4 and 5; the metric prices both at one circuit-wide error, so only the
        # total matters here.
        circuit = km.Circuit("Z_POW q0 0.125\nZ_POW q1 0.0625")
        assert Resources.from_circuit(circuit).num_fine_rotations == 2

    def test_counts_a_qubit_touched_only_by_a_rotation(self, mixed):
        assert mixed.num_qubits == 6

    def test_survives_a_file_round_trip(self, mixed, mixed_circuit, tmp_path):
        path = tmp_path / "c.kmx"
        mixed_circuit.to_file(str(path))
        assert Resources.from_file(path) == mixed

    def test_is_hashable_so_candidates_can_be_deduped(self, mixed, mixed_circuit):
        # frozen=True has to mean it: a search loop dedupes thousands of these in a set.
        assert len({mixed, Resources.from_circuit(mixed_circuit)}) == 1
