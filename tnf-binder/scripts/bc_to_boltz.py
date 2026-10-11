"""Emit Boltz specs for BindCraft designs so they can be scored on the same yardstick.

BindCraft reports AF2 i_pTM on designs it optimised against AF2, which is self-confirming: the
one design checked independently scored i_pTM 0.92 under AF2 but 0.21-0.26 under Boltz. The
campaign designs were all measured under Boltz, so these must be too before any comparison.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.specs import cofold_spec
from tnf_binder.tables import number, read_rows, write_rows
from tnf_binder.target import HUMAN

app = typer.Typer(add_completion=False)
console = Console()

MANIFEST_COLUMNS = ["spec", "design", "kind", "af2_iptm", "af2_ipae", "shape_comp",
                    "n_interface", "rosetta_energy", "positions", "pkas", "sequence"]


@app.command()
def main(
    stats_csv: Path = typer.Argument(..., help="BindCraft's trajectory stats CSV"),
    out_dir: Path = typer.Option(..., "--out-dir", help="where the Boltz specs go"),
    manifest: Path = typer.Option(..., "--manifest", help="CSV mapping spec -> source design"),
    min_iptm: float = typer.Option(0.6, "--min-iptm",
                                   help="AF2 i_pTM floor, applied only to choose what is "
                                        "worth re-folding -- never reported as a result"),
) -> None:
    """Write one Boltz co-folding spec per BindCraft design above the AF2 i_pTM floor."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(stats_csv)
    kept = []
    for row in rows:
        sequence = (row.get("Sequence") or "").strip().upper()
        af2_iptm = number(row, "i_pTM")
        if not sequence or af2_iptm < min_iptm:
            continue
        name = f"BC__{row['Design']}"
        (out_dir / f"{name}.yaml").write_text(cofold_spec(HUMAN, binder=sequence))
        kept.append({"spec": name, "design": row["Design"], "kind": "bindcraft",
                     "af2_iptm": af2_iptm, "af2_ipae": row.get("i_pAE"),
                     "shape_comp": row.get("ShapeComplementarity"),
                     "n_interface": row.get("n_InterfaceResidues"),
                     "rosetta_energy": row.get("Binder_Energy_Score"),
                     "positions": "", "pkas": "", "sequence": sequence})

    write_rows(manifest, kept, MANIFEST_COLUMNS)
    console.print(f"{len(rows)} trajectories, [bold]{len(kept)}[/] with AF2 i_pTM >= {min_iptm}")
    for record in kept:
        console.print(f"  {record['spec']:30s} af2_iptm={record['af2_iptm']} "
                      f"len={len(record['sequence'])}")
    console.print(f"wrote specs -> [cyan]{out_dir}[/], manifest [cyan]{manifest}[/]")


if __name__ == "__main__":
    app()
