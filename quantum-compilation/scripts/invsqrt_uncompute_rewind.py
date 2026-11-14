"""How far an invsqrt circuit's uncompute has rewound its forward pass, and at what Toffoli cost.

Every qubit's value across shots is fingerprinted after every step (one op, or one whole
condition block) from the lookup to the end of the uncompute. A value's birth is the first forward
step holding it; at each uncompute step the latest-born forward value still held marks how far
back the uncompute has got. Matching values rather than whole states is robust to the uncompute
undoing registers on different qubits or in a different order. Findings:
docs/invsqrt-diagnosis.md, "Where the uncompute goes".

    uv run python scripts/invsqrt_uncompute_rewind.py <invsqrt.kmb> --verifiers <dir>
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Annotated

import kickmix as km
import numpy as np
import typer
from rich.progress import track

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.invsqrt import STAGES
from quantum_compilation.submission import Op
from quantum_compilation.verify import load_verifier


def fingerprints(state: np.ndarray) -> np.ndarray:
    """One 64-bit hash per qubit of its bits across shots."""
    return np.array(
        [int.from_bytes(hashlib.blake2b(row.tobytes(), digest_size=8).digest(), "little") for row in np.packbits(state, axis=1)],
        dtype=np.uint64,
    )


def rewind_points(values: np.ndarray, constants: set[int], forward_steps: int) -> np.ndarray:
    """For each step from `forward_steps` on, the latest forward step (> 0) that first produced a
    value still held, or 0 if none is."""
    birth: dict[int, int] = {}
    for s in range(forward_steps + 1):
        for v in values[s].tolist():
            if v not in constants:
                birth.setdefault(v, s)
    return np.array(
        [max((birth.get(v, 0) for v in values[s].tolist()), default=0) for s in range(forward_steps, len(values))]
    )


def simulate(
    circuit: km.Circuit, ops: list[Op], steps: list[tuple[int, int]], x_bits: np.ndarray, x_reg: list[km.q]
) -> tuple[np.ndarray, np.ndarray]:
    """Fired Toffolis per shot in each step, and every qubit's fingerprint before and after each."""
    shots = x_bits.shape[1]
    sim = km.Simulator(batch_size=shots, debug_print=False, seed=0, counting_operations=True)
    sim.write_across_shots(x_reg, x_bits)
    sim.do(Op.circuit(ops[Op.body_start(ops) : steps[0][0]]))
    qubits = [km.q(i) for i in range(circuit.num_qubits)]
    toffolis = np.zeros(len(steps))
    values = [fingerprints(sim.read_across_shots(qubits))]
    for s, (a, b) in enumerate(track(steps, description="steps")):
        sim.clear_operation_statistics()
        sim.do(Op.circuit(ops[a:b]))
        toffolis[s] = sim.operation_statistics().num_toffoli_like / shots
        values.append(fingerprints(sim.read_across_shots(qubits)))
    return toffolis, np.array(values)


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    verifiers: VerifiersOption,
    shots: Annotated[int, typer.Option(help="Uniform inputs, simulated as one batch.")] = 1024,
    report_every: Annotated[float, typer.Option(help="Uncompute Toffolis between report rows.")] = 40.0,
    seed: Annotated[int, typer.Option()] = 0,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    ops = Op.parse(circuit)
    steps = Op.steps(ops, Op.cut_after(ops, STAGES["lookup"]), Op.cut_after(ops, STAGES["clifford"]))
    xs = np.random.default_rng(seed).integers(1, spec.max_x + 1, shots).astype(np.uint64)
    toffolis, values = simulate(circuit, ops, steps, vf.to_bits(xs, spec.width), x_reg)

    cuts = [a for a, _ in steps] + [steps[-1][1]]  # cuts[s]: the op before which step s's state is read
    u0 = cuts.index(Op.cut_after(ops, STAGES["uncompute"]))
    constants = set(fingerprints(np.array([np.zeros(shots, bool), np.ones(shots, bool)])).tolist())
    points = rewind_points(values, constants, u0)
    spent = np.r_[0, np.cumsum(toffolis)]
    typer.echo(f"{'uncompute op':>12} {'uncompute T':>11} | {'rewound to op':>13} {'forward T undone':>16}")
    next_row = 0.0
    for s, p in enumerate(points.tolist(), start=u0):
        cost = spent[s] - spent[u0]
        if cost >= next_row or s == len(values) - 1:
            typer.echo(f"{cuts[s]:>12} {cost:11.1f} | {cuts[p]:>13} {spent[u0] - spent[p]:16.1f}")
            next_row = cost + report_every
    rewound = (cuts[points.min()], cuts[points[0]])  # the forward span the uncompute walks back over
    report_mirror(ops, steps[u0:], toffolis[u0:], rewound)


def report_mirror(
    ops: list[Op], steps: list[tuple[int, int]], toffolis: np.ndarray, rewound: tuple[int, int]
) -> None:
    """A reversal recomputes each AND the forward measured away, so its mirror cost is the forward
    HMR count over the span it rewinds; conditioned Toffolis are measurement-dependent fixups on top."""
    conditioned = sum(t for (a, b), t in zip(steps, toffolis) if any(op.conditioned for op in ops[a:b]))
    measured = sum(op.name == "HMR" for op in ops[slice(*rewound)])
    typer.echo(
        f"uncompute: {toffolis.sum() - conditioned:.1f} unconditioned + {conditioned:.1f} conditioned "
        f"Toffolis per shot; rewound forward ops {rewound[0]}..{rewound[1]} measured away {measured} values"
    )


if __name__ == "__main__":
    typer.run(main)
