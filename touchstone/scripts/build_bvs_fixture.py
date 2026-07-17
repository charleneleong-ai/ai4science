"""Regenerate tests/fixtures/metalpdb_bvs_sites.json — real metalloprotein sites (donor element +
distance) that the recalibrated bond-valence params are checked against (TestBvsCalibration).

Reuses build_bond_valence_params' own fetch + donor_bonds, so the fixture cannot drift from the
recalibration it validates: same shell window, MIN_CN, solvent exclusion.

    uv run python scripts/build_bvs_fixture.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import typer

SCRIPTS = Path(__file__).parent
OUT = SCRIPTS.parent / "tests" / "fixtures" / "metalpdb_bvs_sites.json"
METALS = {"Ni2+": "Ni", "Cu2+": "Cu", "Co2+": "Co"}

spec = importlib.util.spec_from_file_location("bbvp", SCRIPTS / "build_bond_valence_params.py")
bbvp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bbvp)


def main(per_metal: int = 120) -> None:
    out = {}
    for label, symbol in METALS.items():
        sites = bbvp.donor_bonds(bbvp.fetch(symbol), symbol)
        step = max(1, len(sites) // per_metal)  # evenly spaced, not the head (PDB IDs sort by age)
        out[label] = [[[e, round(d, 3)] for e, d in s] for s in sites[::step][:per_metal]]
        print(f"{label}: {len(out[label])} of {len(sites)} real sites", flush=True)
    OUT.write_text(json.dumps(out) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
