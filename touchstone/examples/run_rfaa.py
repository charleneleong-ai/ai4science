"""Rank real designed nickel structures with the geometry verifier.

The adapter is source-agnostic: point it at a directory of RFdiffusionAA backbones or
of LigandMPNN-packed full-atom designs — it parses the metal site either way.

    python examples/run_rfaa.py <design_pdb_dir>

Judged against the real PDB reference. Bare RFdiffusionAA backbones DEFER — no sidechains,
so no real coordination. After LigandMPNN packs His/Asp/Cys sidechains the designs gain
real donors and the scores climb; full trust still needs clean geometry. That progression
is the point — the verifier tracks how real the site is.
"""

import typer
from rich.console import Console
from rich.table import Table

from touchstone import GeometryVerifier, PDBReference, RFdiffusionAdapter, rank

console = Console()
_STYLE = {"trust": "green", "weak": "yellow", "defer": "red"}


def main(design_dir: str) -> None:
    verifier = GeometryVerifier(PDBReference())
    designs = RFdiffusionAdapter(design_dir, pdb_element="NI", metal_label="Ni2+").design("Ni2+")
    table = Table(title=f"touchstone — {len(designs)} nickel designs")
    for col in ("design", "CN", "score", "verdict", "reason"):
        table.add_column(col, justify="left" if col in ("design", "reason") else "right")
    for d, v in rank(designs, verifier):
        table.add_row(d.sequence, str(d.site.coordination_number), f"{v.score:.2e}", v.label, v.reason,
                      style=_STYLE.get(v.label))
    console.print(table)


if __name__ == "__main__":
    typer.run(main)
