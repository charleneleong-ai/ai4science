"""What 6c6a734d's lookup is: its address, its Toffolis, where its CCX gates write, and how
nonlinear its outputs are in the address.

One input per table entry: the mantissa bits above LOW_BITS at an even and an odd exponent
(e = 20, 19, both in range). The lookup spans STAGES["lookup"] to STAGES["lookup cleanup"]
(docs/invsqrt-diagnosis.md).

    uv run python scripts/invsqrt_lookup_decode.py <invsqrt.kmb> --verifiers <dir>
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Annotated

import kickmix as km
import numpy as np
import typer

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.invsqrt import (
    LOOKUP_ADDRESS,
    LOW_BITS,
    MANTISSA_BITS,
    STAGES,
    algebraic_degrees,
    anf,
)
from quantum_compilation.submission import Op
from quantum_compilation.verify import load_verifier

PROBES = 8  # address-register snapshots inside the lookup


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    verifiers: VerifiersOption,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    spec = vf.ProblemSpecInvSqrt.for_n_p(10)
    x_reg, _ = vf.validate_registers(circuit, spec)
    ops = Op.parse(circuit)
    top = (1 << (MANTISSA_BITS - 1 - LOW_BITS)) + np.arange(1 << (MANTISSA_BITS - 1 - LOW_BITS))
    xs = np.r_[top << (LOW_BITS - 1), top << (LOW_BITS - 2)].astype(np.uint64)
    sim = km.Simulator(batch_size=len(xs), debug_print=False, seed=0, counting_operations=True)
    sim.write_across_shots(x_reg, vf.to_bits(xs, spec.width))
    start, stop = Op.cut_after(ops, STAGES["lookup"]), Op.cut_after(ops, STAGES["lookup cleanup"])
    lookup = ops[start:stop]
    everything = [km.q(i) for i in range(circuit.num_qubits)]
    address = [km.q(i) for i in LOOKUP_ADDRESS]

    sim.do(Op.circuit(ops[Op.body_start(ops) : start]))
    shifter = sim.operation_statistics().num_toffoli_like / len(xs)
    typer.echo(f"before the lookup (shifter): {shifter:.2f} Toffolis")
    pre = sim.read_across_shots(everything).astype(np.uint8)
    bits = pre[list(LOOKUP_ADDRESS)]
    order = np.argsort((bits.astype(int) << np.arange(len(bits))[:, None]).sum(axis=0))
    typer.echo(
        f"address q{LOOKUP_ADDRESS[0]}-q{LOOKUP_ADDRESS[-1]} at {start}: "
        f"{np.unique(bits, axis=1).shape[1]} distinct values, degrees {algebraic_degrees(bits[:, order].T)}"
    )

    sim.clear_operation_statistics()
    peak, at = 0, start
    probes = np.linspace(start, stop, PROBES + 1)[1:].astype(int)
    for cut in sorted({Op.cut_after(ops, c) for c in probes}):
        sim.do(Op.circuit(ops[at:cut]))
        at = cut
        peak = max(peak, *algebraic_degrees(sim.read_across_shots(address).astype(np.uint8)[:, order].T))
    toffolis = sim.operation_statistics().num_toffoli_like / len(xs)
    post = sim.read_across_shots(everything).astype(np.uint8)
    outputs = np.flatnonzero((pre != post).any(axis=1))
    names = Counter(op.name for op in lookup)
    conditioned = sum(op.conditioned for op in lookup if op.name == "CCX")
    typer.echo(
        f"lookup: ops {start}-{stop}, {toffolis:.2f} Toffolis per shot ({conditioned} conditioned CCX), "
        f"generic init_lookup {len(xs) - 2}; address degree up to {peak} inside, restored: "
        f"{(post[list(LOOKUP_ADDRESS)] == bits).all()}"
    )

    truth = post[outputs][:, order].T
    nonlinear = (anf(truth).any(axis=1) & (np.bitwise_count(np.arange(len(truth))) >= 2)).sum()
    typer.echo(
        f"{len(outputs)} outputs; degrees {dict(sorted(Counter(algebraic_degrees(truth)).items()))}; "
        f"{nonlinear} of {len(truth) - len(LOOKUP_ADDRESS) - 1} nonlinear monomials used"
    )
    targets = Counter(
        "address" if q in LOOKUP_ADDRESS else "output" if q in outputs else "other"
        for q in (op.qubits[-1] for op in lookup if op.name == "CCX")
    )
    typer.echo(f"{names['CCX']} CCX, targets {dict(targets)}; HMR {names['HMR']}")


if __name__ == "__main__":
    typer.run(main)
