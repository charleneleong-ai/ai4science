"""What the Toffolis inside an invsqrt uncompute's condition blocks are spent on.

Each Toffoli-bearing condition block is costed by simulation and grouped by how many measurement
outcomes its condition bit depends on; blocks keyed on several are listed one by one. Findings:
docs/invsqrt-diagnosis.md, "Where the uncompute goes".

    uv run python scripts/invsqrt_fixups.py <invsqrt.kmb> --verifiers <dir>
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Annotated

import kickmix as km
import numpy as np
import typer

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.invsqrt import STAGES
from quantum_compilation.submission import Op, Resources
from quantum_compilation.verify import load_verifier

def fired_toffolis(
    ops: list[Op], blocks: list[tuple[int, int]], x_bits: np.ndarray, x_reg: list[km.q]
) -> list[float]:
    """Toffolis per shot fired inside each block, simulating everything in order."""
    shots = x_bits.shape[1]
    sim = km.Simulator(batch_size=shots, debug_print=False, seed=0, counting_operations=True)
    sim.write_across_shots(x_reg, x_bits)
    at = Op.body_start(ops)
    fired = []
    for a, b in blocks:
        sim.do(Op.circuit(ops[at:a]))
        sim.clear_operation_statistics()
        sim.do(Op.circuit(ops[a:b]))
        fired.append(sim.operation_statistics().num_toffoli_like / shots)
        at = b
    return fired


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    verifiers: VerifiersOption,
    shots: Annotated[int, typer.Option(help="Uniform inputs, simulated as one batch.")] = 4096,
    seed: Annotated[int, typer.Option()] = 0,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    ops = Op.parse(circuit)
    stop = Op.cut_after(ops, STAGES["clifford"])
    static = {
        (a, b): Resources.from_circuit(Op.circuit(ops[a:b])).num_toffoli_like
        for a, b in Op.steps(ops, Op.cut_after(ops, STAGES["uncompute"]), stop)
        if ops[a].name == "PUSH_CONDITION"
    }
    blocks = [block for block, n in static.items() if n]
    xs = np.random.default_rng(seed).integers(1, spec.max_x + 1, shots).astype(np.uint64)
    fired = fired_toffolis(ops, blocks, vf.to_bits(xs, spec.width), x_reg)

    sources = Op.condition_sources(ops, stop)
    count, static_by, fired_by = Counter(), Counter(), Counter()
    typer.echo("blocks keyed on several outcomes:")
    for (a, b), f in zip(blocks, fired):
        outcomes = len(sources[a])
        kind = "one outcome" if outcomes == 1 else "several"
        count[kind] += 1
        static_by[kind] += static[a, b]
        fired_by[kind] += f
        if outcomes > 1:
            nested = sum(op.name == "PUSH_CONDITION" for op in ops[a + 1 : b])
            typer.echo(
                f"  op {a:>6} {ops[a].condition_bit:>7}: {outcomes:2d} outcomes, {static[a, b]:2d} static, "
                f"{f:4.2f} fired, {nested} nested blocks"
            )
    typer.echo(f"{len(blocks)} condition blocks hold Toffolis; {sum(fired):.1f} fire per shot")
    for kind, f in fired_by.most_common():
        typer.echo(f"  {kind:>11}: {count[kind]:3d} blocks, {static_by[kind]:3d} static, {f:5.1f} fired")


if __name__ == "__main__":
    typer.run(main)
