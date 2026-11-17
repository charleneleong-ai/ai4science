"""In-place lookup synthesis: verified circuits, and free outputs that cost no Toffolis."""

import numpy as np
import pytest
from quantum_compilation.lookup_synth import InPlaceLookup, independent_outputs

INDEX = np.arange(16)
A0, A1 = INDEX & 1, (INDEX >> 1) & 1


@pytest.mark.parametrize(("address_bits", "outputs", "seed"), [(4, 3, 0), (5, 6, 1), (6, 10, 2)])
def test_synthesised_circuit_computes_a_random_table_in_place(address_bits, outputs, seed):
    table = np.random.default_rng(seed).integers(0, 2, (2**address_bits, outputs)).astype(np.uint8)
    synth = InPlaceLookup(table, seed)
    assert synth.run() > 0
    assert synth.verify()


def test_dependent_outputs_cost_no_toffolis():
    nonlinear = A0 & A1
    table = np.stack([nonlinear, np.ones(16, int), A0 ^ 1, nonlinear ^ A1, nonlinear], axis=1)
    assert independent_outputs(table) == [0]
    alone, full = InPlaceLookup(table[:, :1]), InPlaceLookup(table)
    assert alone.run() == full.run()
    assert full.verify()
