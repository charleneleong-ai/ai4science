"""How many Toffolis holding recomputed products on spare clean wires could save in a dqi circuit.

    uv run python -I scripts/dqi_reach.py <circuit.kmb> <work-dir>

The text dump and its `.npz` column cache go in `work-dir` (`CircuitColumns.text_dump`).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import typer

from quantum_compilation.dqi_ir import CircuitColumns


def main(circuit_path: Path, work_dir: Path) -> None:
    columns = CircuitColumns.load(CircuitColumns.text_dump(circuit_path, work_dir))
    found = columns.recomputations()
    held = columns.hold_assignment(found)
    typer.echo(f"{len(found)} recomputed products, {len(held)} holdable on a spare wire")
    typer.echo(f"wires: {dict(Counter(held.values()))}; gaps: {dict(Counter(r.gap for r in held))}")


if __name__ == "__main__":
    typer.run(main)
