"""How many Toffolis holding recomputed products on spare clean wires could save in a dqi circuit.

    uv run python -I scripts/dqi_reach.py <circuit.kmb> <work-dir>

The 1 GB text dump (`<stem>.kmx`) and its `.npz` column cache go in `work-dir`, never beside the
circuit, so a submission folder is left as it is.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import kickmix as km
import typer

from quantum_compilation.dqi_ir import CircuitColumns


def main(circuit_path: Path, work_dir: Path) -> None:
    text_path = work_dir / f"{circuit_path.stem}.kmx"
    if not text_path.exists():
        work_dir.mkdir(parents=True, exist_ok=True)
        km.Circuit.from_file(str(circuit_path)).to_file(str(text_path))
    columns = CircuitColumns.load(text_path)
    found = columns.recomputations()
    held = columns.hold_assignment(found)
    typer.echo(f"{len(found)} recomputed products, {len(held)} holdable on a spare wire")
    typer.echo(f"wires: {dict(Counter(held.values()))}; gaps: {dict(Counter(r.gap for r in held))}")


if __name__ == "__main__":
    typer.run(main)
