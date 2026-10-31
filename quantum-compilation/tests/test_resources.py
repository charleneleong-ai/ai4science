import kickmix as km
import pytest
from quantum_compilation.resources import Resources

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


@pytest.fixture(scope="module")
def mixed_circuit() -> km.Circuit:
    return km.Circuit(MIXED)


@pytest.fixture(scope="module")
def mixed(mixed_circuit) -> Resources:
    return Resources.from_circuit(mixed_circuit)


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
        assert Resources.from_circuit(
            km.Circuit("Z_POW q0 0.125\nZ_POW q1 0.0625")
        ).num_fine_rotations == 2

    def test_counts_a_qubit_touched_only_by_a_rotation(self, mixed):
        assert mixed.num_qubits == 6

    def test_survives_a_file_round_trip(self, mixed, mixed_circuit, tmp_path):
        path = tmp_path / "c.kmx"
        mixed_circuit.to_file(str(path))
        assert Resources.from_file(path) == mixed

    def test_is_hashable_so_candidates_can_be_deduped(self, mixed, mixed_circuit):
        # frozen=True has to mean it: a search loop dedupes thousands of these in a set.
        assert len({mixed, Resources.from_circuit(mixed_circuit)}) == 1
