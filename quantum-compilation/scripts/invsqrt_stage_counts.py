"""Mean Toffolis per shot in each stage of an invsqrt circuit, conditioned gates counted as fired.

Inputs are uniform, as the verifier weights them (its full sweep is too long to rerun per
stage), so the stages sum to its Toffoli count within sampling noise. Stage boundaries are those
decoded for 6c6a734d (docs/invsqrt-diagnosis.md), each moved forward to the next instruction
outside a condition block.

    uv run python scripts/invsqrt_stage_counts.py <invsqrt.kmb> --verifiers <dir>
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

STAGES = {  # name: first instruction, as decoded for 6c6a734d
    "shifter": 0,
    "lookup": 1663,
    "lookup cleanup": 7494,
    "pass 1": 8590,
    "arithmetic": 8970,
    "pass 2": 16430,
    "uncompute": 16750,
    "clifford": 25500,
    "tail": 33000,
}


def stage_circuits(ops: list[Op]) -> dict[str, km.Circuit]:
    starts = [Op.cut_after(ops, s) for s in STAGES.values()]
    stops = [*starts[1:], len(ops)]
    return {name: Op.circuit(ops[a:b]) for name, a, b in zip(STAGES, starts, stops)}


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    verifiers: VerifiersOption,
    batches: Annotated[int, typer.Option()] = 16,
    batch_size: Annotated[int, typer.Option()] = 8192,
    seed: Annotated[int, typer.Option()] = 0,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    stages = stage_circuits(Op.parse(circuit))
    totals = dict.fromkeys(stages, 0)
    sim = km.Simulator(batch_size=batch_size, debug_print=False, seed=0, counting_operations=True)
    rng = np.random.default_rng(seed)
    for _ in track(range(batches), description="batches"):
        xs = rng.integers(1, spec.max_x + 1, batch_size).astype(np.uint64)
        sim.clear_for_shot()
        sim.write_across_shots(x_reg, vf.to_bits(xs, spec.width))
        for name, piece in stages.items():
            sim.clear_operation_statistics()
            sim.do(piece)
            totals[name] += sim.operation_statistics().num_toffoli_like
    simulated = batches * batch_size
    for name, total in totals.items():
        typer.echo(f"{name:>15}: {total / simulated:8.2f} Toffolis per shot")
    typer.echo(f"{'total':>15}: {sum(totals.values()) / simulated:8.2f}")


if __name__ == "__main__":
    typer.run(main)
