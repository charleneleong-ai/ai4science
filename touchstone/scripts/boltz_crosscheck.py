"""Cross-check geometry-trusted designs against an independent Boltz-2 co-fold.

Closes the verifier loop: the geometry verifier trusts a design from its coordination
geometry; Boltz-2 then co-folds the same sequence WITH the metal ion, an independent
predictor. Across Boltz's diffusion samples we compare its confidence (ptm / ligand_iptm)
against the touchstone geometry verdict on Boltz's OWN structures.

Agreement = two independent verifiers both like it. The signal is disagreement — most
informative is `boltz-only`: Boltz folds the complex confidently but leaves the metal
loosely coordinated, so the geometry verifier won't trust it. Even the verifier gets
verified.

    uv run python scripts/boltz_crosscheck.py <boltz_predictions_dir>
"""

from __future__ import annotations

import glob
import json
import re
from collections import defaultdict
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from touchstone import BinderDesign, GeometryVerifier, PDBReference, coordination_site_from_pdb

console = Console()
STYLE = {"both ✓": "green", "boltz-only": "yellow", "geom-only": "cyan", "neither": "red"}


def label(boltz_ok: bool, any_trust: bool, any_coord: bool) -> str:
    if boltz_ok and any_trust:
        return "both ✓"
    if boltz_ok:
        return "boltz-only"  # confident fold, but geometry not trustworthy
    if any_trust or any_coord:
        return "geom-only"
    return "neither"


def main(predictions_dir: str, metal_element: str = "NI", metal_label: str = "Ni2+") -> None:
    verifier = GeometryVerifier(PDBReference())
    by_design: dict[str, list] = defaultdict(list)
    for pdb in sorted(glob.glob(f"{predictions_dir}/**/*_model_*.pdb", recursive=True)):
        stem = Path(pdb).stem
        design = re.match(r"(.+?)_model_\d+", stem).group(1)
        conf = json.loads((Path(pdb).parent / f"confidence_{stem}.json").read_text())
        try:
            site = coordination_site_from_pdb(pdb, metal_element, metal_label)
            v = verifier.verify(BinderDesign(design, site, "boltz", float("nan")))
            by_design[design].append((site.coordination_number, v.trust, not v.ood, conf))
        except ValueError:
            by_design[design].append((0, False, False, conf))

    table = Table(title=f"Boltz-2 co-fold vs touchstone geometry ({metal_label}) — over diffusion samples")
    for col in ("design", "n", "best CN", "coord", "trust", "ptm", "lig_iptm", "agree"):
        table.add_column(col, justify="left" if col == "design" else "right")
    for design in sorted(by_design):
        rows = by_design[design]
        any_trust = any(r[1] for r in rows)
        coord = sum(1 for r in rows if r[2])
        ptm = max(r[3]["ptm"] for r in rows)
        lig = max(r[3]["ligand_iptm"] for r in rows)
        agree = label(ptm > 0.7 and lig > 0.6, any_trust, coord > 0)
        table.add_row(design, str(len(rows)), str(max(r[0] for r in rows)), f"{coord}/{len(rows)}",
                      str(any_trust), f"{ptm:.2f}", f"{lig:.2f}", agree, style=STYLE.get(agree))
    console.print(table)


if __name__ == "__main__":
    typer.run(main)
