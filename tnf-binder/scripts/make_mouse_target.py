"""Build a Boltz spec for the mouse TNF-alpha trimer.

Objective 2 wants one binder that recognises human AND mouse TNF-alpha. BindCraft2 supports
multi-target campaigns ("target": [human, mouse]), which optimises both interfaces directly
instead of designing against human and trusting that 79% identity carries the binder across.
That needs a mouse structure, and no mouse TNF-alpha crystal structure is being used here, so
predict the trimer from sequence.

Mouse soluble domain from UniProt P06804, aligned to the human construct boundaries.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.specs import cofold_spec
from tnf_binder.target import MOUSE

app = typer.Typer(add_completion=False)
console = Console()


@app.command()
def main(out_dir: Path = typer.Argument(..., help="directory to write the spec into")) -> None:
    """Emit a three-protomer mouse TNF-alpha spec for Boltz to fold."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "mouse_tnf_trimer.yaml"
    path.write_text(cofold_spec(MOUSE))
    console.print(f"mouse soluble domain [bold]{len(MOUSE)}[/] aa (UniProt P06804)")
    console.print(f"wrote [cyan]{path}[/] (3 protomers)")


if __name__ == "__main__":
    app()
