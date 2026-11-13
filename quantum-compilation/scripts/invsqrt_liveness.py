"""Which qubits of an invsqrt circuit hold a value at each step, and which of those sit parked.

Method and findings: docs/invsqrt-diagnosis.md, "Qubit liveness".

    uv run python scripts/invsqrt_liveness.py <invsqrt.kmb> <out.npz> --verifiers <dir>
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import kickmix as km
import numpy as np
import typer
from rich.progress import track

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.submission import Op
from quantum_compilation.verify import load_verifier


class Steps:
    """The circuit cut into steps: one unconditioned op, or one whole condition block."""

    def __init__(self, circuit: km.Circuit) -> None:
        ops = Op.parse(circuit)
        self.num_qubits = circuit.num_qubits
        self.groups: list[list[Op]] = []
        self.instructions: list[int] = []  # each step's first op, indexed like `str(circuit)`
        current: list[Op] = []
        for i, op in enumerate(ops):
            if op.name in Op.DECLARATIONS:
                continue
            if not current:
                self.instructions.append(i)
            current.append(op)
            if op.depth == 0:
                self.groups.append(current)
                current = []

    def pieces(self) -> list[km.Circuit]:
        return [km.Circuit("\n".join(op.text for op in group)) for group in self.groups]

    def real_uses(self) -> np.ndarray:
        """`used[s, q]`: an unconditioned gate other than HMR uncompute acts on q at step s."""
        used = np.zeros((len(self.groups), self.num_qubits), dtype=bool)
        for s, group in enumerate(self.groups):
            for op in group:
                if not op.conditioned and op.name != "HMR":
                    used[s, list(op.qubits)] = True
        return used


def sample_inputs(max_x: int, count: int, seed: int) -> np.ndarray:
    """Every x up to 4096, both neighbours of each power of two, and the rest uniform.

    Not the verifier's `sample_shots` sweep: occupancy is a max over inputs, so the edge cases
    of the normalising shift matter more than matching its order.
    """
    edges = [v for k in range(1, max_x.bit_length()) for v in (2**k - 1, 2**k, 2**k + 1)]
    fixed = np.unique(np.r_[np.arange(1, 4097), edges, max_x].clip(1, max_x))
    rest = np.random.default_rng(seed).integers(1, max_x + 1, max(count - len(fixed), 0))
    return np.r_[fixed, rest].astype(np.uint64)


def occupancy(
    steps: Steps, x_bits: np.ndarray, x_reg: list[km.q], batch_size: int
) -> np.ndarray:
    """`occupied[s, q]`: qubit q is nonzero, in some shot, before or after step s."""
    pieces = steps.pieces()
    qubits = [km.q(i) for i in range(steps.num_qubits)]
    occupied = np.zeros((len(pieces), len(qubits)), dtype=bool)
    sim = km.Simulator(batch_size=batch_size, debug_print=False, seed=0)
    for start in range(0, x_bits.shape[1], batch_size):
        batch = np.resize(x_bits[:, start : start + batch_size], (len(x_reg), batch_size))
        sim.clear_for_shot()
        sim.write_across_shots(x_reg, batch)  # a short last batch repeats its inputs
        before = sim.read_across_shots(qubits).any(axis=1)
        for s, piece in enumerate(track(pieces, description=f"inputs from {start}")):
            sim.do(piece)
            after = sim.read_across_shots(qubits).any(axis=1)
            occupied[s] |= before | after
            before = after
    return occupied


def runs(mask: np.ndarray, min_steps: int) -> list[tuple[int, int]]:
    """`[start, stop)` step ranges where `mask` holds for at least `min_steps` steps."""
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(np.int8), 0]))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2]) if b - a >= min_steps]


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    out: Annotated[Path, typer.Argument(help="Where to write the occupancy matrix (.npz).")],
    verifiers: VerifiersOption,
    samples: Annotated[int, typer.Option(help="Inputs to simulate.")] = 16384,
    batch_size: Annotated[int, typer.Option()] = 4096,
    seed: Annotated[int, typer.Option()] = 0,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    xs = sample_inputs(spec.max_x, samples, seed)
    steps = Steps(circuit)
    occupied = occupancy(steps, vf.to_bits(xs, spec.width), x_reg, batch_size)
    used = steps.real_uses()
    instructions = np.array(steps.instructions)
    np.savez_compressed(out, occupied=occupied, used=used, instructions=instructions, xs=xs)

    count = occupied.sum(axis=1)
    peak = int(count.max())
    typer.echo(f"{len(xs)} inputs; {steps.num_qubits} qubits declared; peak occupancy {peak}")
    for a, b in runs(count == peak, min_steps=50):
        parked = np.flatnonzero(occupied[a:b].all(axis=0) & ~used[a:b].any(axis=0))
        typer.echo(
            f"  peak run: instructions {instructions[a]}..{instructions[b - 1]} ({b - a} steps), "
            f"{len(parked)} qubits parked throughout: {' '.join(f'q{q}' for q in parked)}"
        )


if __name__ == "__main__":
    typer.run(main)
