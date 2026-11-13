"""Phase error of an invsqrt circuit on every input x, through the verifier's own shot loop.

Mirrors the verifier's phase check only (not its x-restored or scratch-clean checks), keeping
each x's error where the verifier keeps just the worst ratio and the smallest margin.

    uv run python scripts/invsqrt_error_profile.py <invsqrt.kmb> <out.npz> --verifiers <dir>
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Annotated

import kickmix as km
import numpy as np
import typer
from rich.progress import track

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.cost import (
    PHASE_UNITS_PER_TURN,
    rotation_eps,
    toffoli_equivalents,
)
from quantum_compilation.submission import Resources
from quantum_compilation.verify import load_verifier


def profile(
    circuit: km.Circuit, vf: ModuleType, batch_size: int, shots_per_x: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per x: worst |error|, first-shot signed error, whether shots disagree, and `tol(x)`.

    All in units of 2^-64 turns, kept as running per-x summaries so memory does not grow
    with the number of shots.
    """
    spec = vf.ProblemSpecInvSqrt.for_n_p(10, shots_per_x=shots_per_x)
    x_reg, _ = vf.validate_registers(circuit, spec)
    target, tol = vf.reference_values(
        spec, np.arange(1, spec.max_x + 1, dtype=np.uint64)
    )
    worst = np.zeros(spec.max_x, dtype=np.uint64)
    first = np.zeros(spec.max_x, dtype=np.int64)
    varies = np.zeros(spec.max_x, dtype=bool)
    sim = km.Simulator(batch_size=batch_size, debug_print=False, seed=0)
    for start in track(range(0, spec.num_shots, batch_size), description="shots"):
        # The simulator wants full batches: run one, then drop shots past the end.
        xs = vf.sample_shots(spec, start, start + batch_size)
        sim.clear_for_shot()
        sim.write_across_shots(x_reg, vf.to_bits(xs, spec.width))
        sim.do(circuit)
        xs = xs[: spec.num_shots - start]
        theta = np.array(
            [vf.half_turns_to_units(sim.read_shot_phase(i)) for i in range(len(xs))],
            dtype=np.uint64,
        )
        index = (xs - 1).astype(np.int64)
        error = (theta - target[index]).view(np.int64)
        repeat = np.arange(start, start + len(xs)) >= spec.max_x  # shot after x's first
        first[index[~repeat]] = error[~repeat]
        varies[index[repeat]] |= error[repeat] != first[index[repeat]]
        np.maximum.at(worst, index, np.abs(error).view(np.uint64))
    return worst, first, varies, tol


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    out: Annotated[
        Path, typer.Argument(help="Where to write the per-x summaries (.npz).")
    ],
    verifiers: VerifiersOption,
    shots_per_x: Annotated[
        int,
        typer.Option(
            help="Shots per x; the verifier's 64 would take hours single-process."
        ),
    ] = 1,
    batch_size: Annotated[int, typer.Option()] = 256,
) -> None:
    vf = load_verifier(verifiers, "invsqrt")
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    worst, first, varies, tol = profile(circuit, vf, batch_size, shots_per_x)
    np.savez_compressed(out, worst=worst, first=first, varies=varies, tol=tol)

    n_rot = Resources.from_circuit(circuit).num_fine_rotations
    ratio = worst / tol
    margin = tol.astype(np.int64) - worst.astype(np.int64)
    order = np.argsort(margin)
    top = int(ratio.argmax())
    typer.echo(f"shots varying per x: {int(varies.sum())} of {len(tol)}")
    typer.echo(f"worst ratio {ratio[top]:.5f} at x={top + 1}")
    for k in (0, 1, 3, 10, 30, 100, 1000):
        i = int(order[k])
        eps = rotation_eps(int(margin[i]), n_rot)
        typer.echo(
            f"bind after fixing {k:>4} x: x={i + 1:>8} "
            f"margin {margin[i] / PHASE_UNITS_PER_TURN:.3e} turns (ratio {ratio[i]:.4f}) "
            f"-> eps {eps:.3e}, rotations {toffoli_equivalents(0, 0, n_rot, eps):.1f} C"
        )
    signed = first / tol.astype(np.float64)
    typer.echo(
        f"signed error/tol: mean {signed.mean():+.4f}, sd {signed.std():.4f}; "
        f"ratio > 0.5 for {int((ratio > 0.5).sum())} x, > 0.9 for {int((ratio > 0.9).sum())}"
    )


if __name__ == "__main__":
    typer.run(main)
