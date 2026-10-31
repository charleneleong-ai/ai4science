import kickmix as km
import pytest

from quantum_compilation.cost import (
    PROVISIONAL_CALIBRATION,
    MetricCalibration,
    UnpricedRotationError,
    t_count,
    toffoli_equivalents,
)
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

PRICED = MetricCalibration(
    reference_costs=PROVISIONAL_CALIBRATION.reference_costs,
    rotation_t_cost={4: 10.0, 5: 14.0},
    provenance="test",
    measured=True,
)


@pytest.fixture(scope="module")
def mixed_circuit() -> km.Circuit:
    return km.Circuit(MIXED)


@pytest.fixture(scope="module")
def mixed(mixed_circuit) -> Resources:
    return Resources.from_circuit(mixed_circuit)


class TestResourceExtraction:
    """Width and non-Clifford content, read off the circuit by kickmix."""

    def test_counts_toffoli_likes_including_cswap(self, mixed):
        assert mixed.num_toffoli_like == 2

    def test_separates_t_level_from_finer_rotations(self, mixed):
        # The S gate is Clifford, so it is neither a T nor a fine rotation.
        assert mixed.num_t_like == 4
        assert mixed.fine_rotations == ((4, 1),)
        assert mixed.num_fine_rotations == 1

    def test_counts_a_qubit_touched_only_by_a_rotation(self, mixed):
        assert mixed.num_qubits == 6

    def test_survives_a_file_round_trip(self, mixed, mixed_circuit, tmp_path):
        path = tmp_path / "c.kmx"
        mixed_circuit.to_file(str(path))
        assert Resources.from_file(path) == mixed

    def test_is_hashable_so_candidates_can_be_deduped(self, mixed, mixed_circuit):
        # frozen=True has to mean it: a search loop dedupes thousands of these in a set.
        assert len({mixed, Resources.from_circuit(mixed_circuit)}) == 1


class TestToffoliEquivalents:
    """C = Toffolis + T/4, where the price of a fine rotation comes from the calibration."""

    def test_refuses_to_invent_a_price_for_fine_rotations(self, mixed):
        with pytest.raises(UnpricedRotationError, match=r"level\(s\) \[4\]"):
            toffoli_equivalents(mixed, PROVISIONAL_CALIBRATION)

    def test_prices_fine_rotations_from_the_calibration(self, mixed):
        # 2 Toffolis + (4 T + one level-4 rotation at 10 T) / 4
        assert toffoli_equivalents(mixed, PRICED) == pytest.approx(2 + 14 / 4)

    def test_a_circuit_without_fine_rotations_needs_no_price(self):
        clifford = Resources.from_circuit(km.Circuit("CCX q0 q1 q2\nZ_POW q0 0.5"))
        assert toffoli_equivalents(clifford, PROVISIONAL_CALIBRATION) == pytest.approx(
            1.0
        )

    def test_prices_each_level_separately(self):
        res = Resources.from_circuit(km.Circuit("Z_POW q0 0.125\nZ_POW q1 0.0625"))
        assert res.fine_rotations == ((4, 1), (5, 1))
        assert t_count(res, PRICED) == pytest.approx(24.0)

    def test_names_only_the_levels_the_calibration_is_missing(self):
        res = Resources.from_circuit(km.Circuit("Z_POW q0 0.125\nZ_POW q1 0.03125"))
        partial = MetricCalibration(
            PROVISIONAL_CALIBRATION.reference_costs, {4: 10.0}, "partial", measured=True
        )
        with pytest.raises(UnpricedRotationError, match=r"level\(s\) \[6\]"):
            t_count(res, partial)
