"""Regenerate tests/fixtures/metalpdb_real_sites.json — the real metalloprotein sites the geometry
prior is calibrated against (TestShippedPriorCalibration).

Reuses `build_metalpdb_reference`'s own `fetch` + `donor_shells`, so the fixture cannot drift from
the prior it checks: same shell window, same MIN_CN, same solvent exclusion. If those filters ever
diverge, the calibration number means nothing.

    uv run python scripts/build_real_site_fixture.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import typer

SCRIPTS = Path(__file__).parent
OUT = SCRIPTS.parent / "tests" / "fixtures" / "metalpdb_real_sites.json"
METALS = {"Ni": "Ni2+", "Cu": "Cu2+", "Co": "Co2+"}

spec = importlib.util.spec_from_file_location("bmp", SCRIPTS / "build_metalpdb_reference.py")
bmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bmp)


def main(per_metal: int = 120) -> None:
    out = {}
    for symbol, label in METALS.items():
        shells = [sorted(round(d, 3) for d in s) for s in bmp.donor_shells(bmp.fetch(symbol), symbol)]
        step = max(1, len(shells) // per_metal)  # evenly spaced, not the head (PDB IDs sort by age)
        out[label] = [{"donors": s} for s in shells[::step][:per_metal]]
        print(f"{label}: {len(out[label])} of {len(shells)} real sites", flush=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
