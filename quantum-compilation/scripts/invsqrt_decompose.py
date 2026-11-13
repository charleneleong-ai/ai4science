"""What an invsqrt circuit's lookup reads and what its two kick passes add, from simulated state.

Method and findings: docs/invsqrt-diagnosis.md, "What the circuit computes".

    uv run python scripts/invsqrt_decompose.py <invsqrt.kmb> <out.npz> --verifiers <dir>
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Annotated

import kickmix as km
import numpy as np
import typer

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.cost import PHASE_UNITS_PER_TURN
from quantum_compilation.invsqrt import MANTISSA_BITS, normalised, sample_inputs
from quantum_compilation.submission import Op
from quantum_compilation.verify import load_verifier

LOW_BITS = 14  # mantissa bits below the lookup address, read from 6c6a734d
PARITY_QUBIT = 23  # holds the exponent's parity through the lookup (see the doc)
MIN_FIT_SURPLUS = 3  # samples beyond the coefficient count before a segment is fitted


class Kicks:
    """The fine rotations: where they sit, which qubit carries the bit, and their Z exponent."""

    def __init__(self, ops: list[Op]) -> None:
        fine = [(i, op) for i, op in enumerate(ops) if op.is_fine_rotation]
        self.at = np.array([i for i, _ in fine])
        self.qubits = [op.qubits[0] for _, op in fine]
        self.turns = np.array([float(op.args[1]) for _, op in fine]) / 2  # a kicked 1 adds a/2

    def passes(self) -> list[np.ndarray]:
        """Kick indices split at the widest gap between consecutive kicks: the two passes."""
        cut = int(np.diff(self.at).argmax()) + 1
        return [np.arange(cut), np.arange(cut, len(self.at))]


def simulate(
    circuit: km.Circuit, ops: list[Op], kicks: Kicks, x_bits: np.ndarray, x_reg: list[km.q], lookup_end: int
) -> tuple[np.ndarray, np.ndarray]:
    """All qubits at `lookup_end`, then the bit on the kick qubit just before each fine rotation."""
    if kicks.at[0] < lookup_end or any(ops[i - 1].depth for i in (lookup_end, *kicks.at)):
        raise ValueError("cuts must follow the lookup and fall outside condition blocks")
    sim = km.Simulator(batch_size=x_bits.shape[1], debug_print=False, seed=0)
    sim.write_across_shots(x_reg, x_bits)
    start = sum(op.name in Op.DECLARATIONS for op in ops)

    def run_to(stop: int) -> None:
        nonlocal start
        sim.do(km.Circuit("\n".join(op.text for op in ops[start:stop])))
        start = stop

    run_to(lookup_end)
    state = sim.read_across_shots([km.q(i) for i in range(circuit.num_qubits)]).T
    kicked = []
    for at, qubit in zip(kicks.at, kicks.qubits):
        run_to(int(at))
        kicked.append(sim.read_across_shots(km.q(qubit)))
    return state, np.array(kicked)


def collisions(keys: np.ndarray, outputs: np.ndarray) -> int:
    """Keys seen with more than one output row."""
    _, key_ids = np.unique(keys, axis=0, return_inverse=True)
    _, out_ids = np.unique(outputs, axis=0, return_inverse=True)
    pairs = np.unique(np.c_[key_ids.ravel(), out_ids.ravel()], axis=0)
    return len(pairs) - int(key_ids.max()) - 1


def worst_fit(
    u: np.ndarray, y: np.ndarray, tolerance: np.ndarray, key: np.ndarray, degree: int
) -> float:
    """Worst residual of a per-key polynomial in `u`, in units of each input's tolerance."""
    worst = 0.0
    for k in np.unique(key):
        m = key == k
        if m.sum() > degree + MIN_FIT_SURPLUS:
            fit = np.polyval(np.polyfit(u[m], y[m], degree), u[m])
            worst = max(worst, float((np.abs(fit - y[m]) / tolerance[m]).max()))
    return worst


def report_lookup(state: np.ndarray, xs: np.ndarray, x_reg: list[km.q]) -> None:
    held = set(range(state.shape[1])) - {q.id for q in x_reg} - {PARITY_QUBIT}
    outputs = state[:, [q for q in sorted(held) if state[:, q].any()]]  # erased scratch reads 0
    mantissa, exponent = normalised(xs)
    typer.echo(f"lookup: {outputs.shape[1]} output qubits; key collisions by address width")
    for bits in (6, 7, 8):
        address = mantissa >> (MANTISSA_BITS - 1 - bits) & (2**bits - 1)
        keys = np.c_[address, exponent % 2]
        typer.echo(f"  {bits} mantissa bits + exponent parity: {collisions(keys, outputs)}")


def report_passes(terms: list[np.ndarray], xs: np.ndarray, target: np.ndarray) -> None:
    for x in (100, 1000, 4096):
        j = int(np.flatnonzero(xs == x)[0])
        split = " + ".join(f"{t[j]:.9f}" for t in terms)
        typer.echo(f"x={x}: target {target[j]:.9f} = passes {split}")


def report_fits(terms: list[np.ndarray], xs: np.ndarray, tolerance: np.ndarray) -> None:
    """The last pass, undone of its 2^-floor(e/2) scale, as a polynomial per segment."""
    big = xs >= 2 ** (LOW_BITS + 1)  # every low mantissa bit is an input bit
    mantissa, exponent = normalised(xs[big])
    scale = 2.0 ** (exponent // 2)
    key = (mantissa >> LOW_BITS) | (exponent % 2) << (MANTISSA_BITS - LOW_BITS)
    u = (mantissa & (2**LOW_BITS - 1)) / 2**LOW_BITS
    for degree in (1, 2, 3):
        worst = worst_fit(u, terms[-1][big] * scale, tolerance[big] * scale, key, degree)
        typer.echo(
            f"last pass as a degree-{degree} polynomial in the low {LOW_BITS} bits per segment: "
            f"worst residual {worst:.3f} x tolerance"
        )


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    out: Annotated[Path, typer.Argument(help="Where to write the kicked bits (.npz).")],
    verifiers: VerifiersOption,
    lookup_end: Annotated[
        int, typer.Option(help="Instruction after the lookup (7,494 in 6c6a734d; see the doc).")
    ] = 7494,
    samples: Annotated[int, typer.Option(help="Inputs to simulate.")] = 16384,
    seed: Annotated[int, typer.Option()] = 0,
) -> None:
    vf: ModuleType = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    xs = sample_inputs(spec.max_x, samples, seed)
    ops = Op.parse(circuit)
    kicks = Kicks(ops)
    state, kicked = simulate(circuit, ops, kicks, vf.to_bits(xs, spec.width), x_reg, lookup_end)
    np.savez_compressed(out, xs=xs, state=state, kicked=kicked, kick_at=kicks.at)

    target, tol = vf.reference_values(spec, xs)
    terms = [(kicked[p] * kicks.turns[p, None]).sum(axis=0) for p in kicks.passes()]
    report_lookup(state, xs, x_reg)
    report_passes(terms, xs, target / PHASE_UNITS_PER_TURN)
    report_fits(terms, xs, tol / PHASE_UNITS_PER_TURN)


if __name__ == "__main__":
    typer.run(main)
