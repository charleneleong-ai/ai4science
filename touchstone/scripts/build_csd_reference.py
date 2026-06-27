"""Build touchstone's CSD reference geometry from the Cambridge Structural Database.

Mirrors build_pdb_reference.py but draws on the CSD's small-molecule metal–organic
coordination — sharper donor-geometry priors for chelator-style designs (the
metal-recovery use case) — instead of the PDB's protein sites. Writes the empirical
mean/std metal–donor bond length + modal coordination number + observed CN range to
src/touchstone/data/csd_reference.json, the file CSDReference loads. Same schema as
pdb_reference.json, so the verifier is unchanged.

Requires the CSD Python API (`ccdc`) and a valid CSD licence (gated — this is why the
data file is not committed):

    conda install -c "https://conda.ccdc.cam.ac.uk" csd-python-api   # licence required
    uv run python scripts/build_csd_reference.py --max-hits 2000
"""

from __future__ import annotations

import json
import statistics
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import typer

OUT = Path(__file__).parent.parent / "src" / "touchstone" / "data" / "csd_reference.json"
METALS = {"Ni": "Ni2+", "Cu": "Cu2+", "Co": "Co2+"}  # CSD element symbol → verifier label
DONORS = {"N", "O", "S"}
SHELL = (1.0, 2.8)  # first-shell distance window, Å
MIN_CN = 3  # ignore adventitious ions with < 3 contacts


def _entries_with(metal_symbol: str, max_hits: int) -> list[str]:
    """CSD refcodes whose structures contain the metal (element substructure search)."""
    from ccdc.search import SMARTSSubstructure, SubstructureSearch

    search = SubstructureSearch()
    search.add_substructure(SMARTSSubstructure(f"[{metal_symbol}]"))
    return [hit.identifier for hit in search.search(max_hit_structures=max_hits)]


def _sites(molecule, metal_symbol: str):
    """Per metal atom in a CSD molecule, the first-shell donor bond lengths (Å)."""
    atoms = [(a.atomic_symbol, np.array(a.coordinates)) for a in molecule.atoms if a.coordinates]
    metals = [xyz for el, xyz in atoms if el == metal_symbol]
    donors = [xyz for el, xyz in atoms if el in DONORS]
    for m in metals:
        d = np.linalg.norm(np.array(donors) - m, axis=1) if donors else np.empty(0)
        shell = d[(d > SHELL[0]) & (d <= SHELL[1])]
        if len(shell) >= MIN_CN:
            yield shell


def build_metal(metal_symbol: str, label: str, max_hits: int) -> dict:
    from ccdc import io

    bonds: list[float] = []
    counts: list[int] = []
    refcodes = _entries_with(metal_symbol, max_hits)
    with io.EntryReader("CSD") as reader:
        for refcode in refcodes:
            try:
                mol = reader.entry(refcode).molecule
            except Exception:
                continue
            for shell in _sites(mol, metal_symbol):
                bonds.extend(shell.tolist())
                counts.append(len(shell))
    return {
        "metal": label,
        "coordination_number": Counter(counts).most_common(1)[0][0],
        "cn_range": [int(np.percentile(counts, 10)), int(np.percentile(counts, 90))],
        "bond_length_mean": round(statistics.fmean(bonds), 3),
        "bond_length_std": round(statistics.pstdev(bonds), 3),
        "source": f"CSD ≤{max_hits} hits, {len(counts)} sites / {len(bonds)} bonds, pulled {date.today()}",
    }


def main(max_hits: int = 2000) -> None:
    table = {}
    for metal_symbol, label in METALS.items():
        print(f"building {label} from CSD element {metal_symbol} ...", flush=True)
        table[label] = build_metal(metal_symbol, label, max_hits)
        print(f"  {table[label]}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
