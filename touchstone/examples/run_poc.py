"""POC: generator → geometry verifier → ranked trust/defer.

Runs MockGenerator + MockReference end-to-end and shows the extreme-condition
(acidic-leachate) angle. Swap MockGenerator → RFdiffusionAdapter and MockReference →
PDBReference / CSD-Mogul and the loop below does not change.

    python examples/run_poc.py
"""

import typer
from rich.console import Console
from rich.table import Table

from touchstone import GeometryVerifier, MockGenerator, design_and_rank, under_leachate

TARGET = "Ni2+"
console = Console()
_STYLE = {"trust": "green", "weak": "yellow", "defer": "red"}


def main() -> None:
    generator = MockGenerator(seed=0)
    verifier = GeometryVerifier()  # MockReference by default
    ranked = design_and_rank(generator, verifier, TARGET, n=5)

    table = Table(title=f"touchstone POC — {TARGET}")
    for col in ("design", "score", "verdict", "reason"):
        table.add_column(col, justify="left" if col in ("design", "reason") else "right")
    for design, verdict in ranked:
        table.add_row(design.sequence, f"{verdict.score:.3f}", verdict.label, verdict.reason,
                      style=_STYLE.get(verdict.label))
    console.print(table)

    # The extreme-condition angle: the best design, re-judged under acidic leachate.
    best, _ = ranked[0]
    stressed = verifier.verify(under_leachate(best, bond_stretch=0.6))
    console.print(f"\n[bold]{best.sequence}[/] under acidic leachate → "
                  f"score {stressed.score:.3f}, [red]ood={stressed.ood}[/] ({stressed.reason})")


if __name__ == "__main__":
    typer.run(main)
