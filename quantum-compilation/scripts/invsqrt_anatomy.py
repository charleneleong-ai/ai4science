"""Cost anatomy of one invsqrt circuit: where its Toffolis sit and what its rotations cost.

Reads the circuit statically (no simulation) and, optionally, the verifier's own row for it, so
the rotation price uses the synthesis error the verifier actually chose.

    uv run python scripts/invsqrt_anatomy.py <invsqrt.kmb> [--rows rows.jsonl] [--window 1500]
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Annotated

import kickmix as km
import typer

from quantum_compilation.cost import T_PER_TOFFOLI_EQUIVALENT, rotation_t_count
from quantum_compilation.submission import Resources
from quantum_compilation.verify import VerificationRun

TOFFOLI_LIKE = ("CCX", "CCZ", "CSWAP")
TRACKED = (*TOFFOLI_LIKE, "HMR", "Z_POW")


def timeline(ops: list[list[str]], window: int) -> list[tuple[int, Counter[str], int]]:
    rows, current, conditioned, depth = [], Counter(), 0, 0
    for i, op in enumerate(ops, start=1):
        name = op[0]
        depth += (name == "PUSH_CONDITION") - (name == "POP_CONDITION")
        if name in TRACKED:
            current[name] += 1
            conditioned += bool(depth) and name in TOFFOLI_LIKE
        if i % window == 0 or i == len(ops):
            rows.append((i, current, conditioned))
            current, conditioned = Counter(), 0
    return rows


def fine_rotations(ops: list[list[str]]) -> list[tuple[str, float]]:
    """`Z_POW`s whose exponent is not a multiple of T's (Z^0.25), i.e. the ones needing synthesis."""
    return [
        (op[1], a)
        for op in ops
        if op[0] == "Z_POW" and not (4 * (a := float(op[2]))).is_integer()
    ]


def main(
    circuit_path: Annotated[Path, typer.Argument(help="The invsqrt .kmb file.")],
    rows: Annotated[
        Path | None, typer.Option(help="rows.jsonl from `vqcc verify`.")
    ] = None,
    window: Annotated[int, typer.Option(help="Instructions per timeline row.")] = 1500,
) -> None:
    circuit = km.Circuit.from_file(str(circuit_path), format="kmb")
    ops = [line.split() for line in str(circuit).splitlines() if line.strip()]
    typer.echo(
        f"{'upto':>6} " + " ".join(f"{k:>6}" for k in TRACKED) + "  conditioned_toffoli"
    )
    for upto, counts, conditioned in timeline(ops, window):
        typer.echo(
            f"{upto:>6} "
            + " ".join(f"{counts[k]:>6}" for k in TRACKED)
            + f"  {conditioned:>19}"
        )

    angles = fine_rotations(ops)
    expected = Resources.from_circuit(circuit).num_fine_rotations
    if len(angles) != expected:
        raise typer.BadParameter(
            f"found {len(angles)} fine rotations, kickmix counts {expected}"
        )
    exponents = Counter(round(math.log2(a)) for _, a in angles)
    typer.echo(
        f"\n{len(angles)} fine rotations on qubits {sorted({q for q, _ in angles})}; "
        f"exponents {min(exponents)}..{max(exponents)}, "
        f"repeated: {sorted(e for e, n in exponents.items() if n > 1)}"
    )

    if rows is None:
        return
    detail = VerificationRun.read_rows(rows)["invsqrt"]["detail"]
    t_per_rotation = rotation_t_count(detail["rotation_eps"])
    c = detail["toffoli_equivalent"]
    rotation_c = detail["rotation_count"] * t_per_rotation / T_PER_TOFFOLI_EQUIVALENT
    typer.echo(
        f"verifier: Q {detail['num_qubits']}, C {c:.2f} = Toffoli {detail['toffoli_count']:.2f} "
        f"+ rotations {rotation_c:.2f} ({rotation_c / c:.1%}); {t_per_rotation:.2f} T per rotation at "
        f"eps {detail['rotation_eps']:.3g}; worst error/tolerance {detail['worst_error_over_tolerance']:.5f} "
        f"at x={detail['worst_x']}"
    )


if __name__ == "__main__":
    typer.run(main)
