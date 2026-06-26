"""Build touchstone's reference metal-coordination geometry from the public PDB.

Pulls high-resolution structures containing each metal from RCSB, measures every
first-shell metal–donor (N/O/S) bond, and writes the empirical mean/std + modal
coordination number + observed CN range to src/touchstone/data/pdb_reference.json —
which PDBReference loads. Reproducible: re-run to refresh.

    uv run python scripts/build_pdb_reference.py
    uv run python scripts/build_pdb_reference.py --n-structures 60 --max-res 1.5
"""

from __future__ import annotations

import json
import statistics
import urllib.parse
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import typer

from touchstone.geometry.parse import DONOR_ELEMENTS

OUT = Path(__file__).parent.parent / "src" / "touchstone" / "data" / "pdb_reference.json"
METALS = {"NI": "Ni2+", "CU": "Cu2+"}  # PDB comp_id → verifier label
SHELL = (1.0, 2.8)  # first-shell distance window
MIN_CN = 3  # ignore adventitious ions with < 3 contacts


def _search(comp_id: str, n_structures: int, max_res: float) -> list[str]:
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text_chem", "parameters": {
                    "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                    "operator": "exact_match", "value": comp_id}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less", "value": max_res}},
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": n_structures},
                            "results_content_type": ["experimental"]},
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(query))
    res = json.load(urllib.request.urlopen(url, timeout=30))
    return [hit["identifier"] for hit in res["result_set"]]


def _atoms(pdb_text: str):
    """Yield (element, xyz, altloc) for ATOM/HETATM records."""
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[76:78].strip().upper() or line[12:16].strip().lstrip("0123456789")[:1].upper()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        yield el, xyz, line[16]


def _sites(pdb_text: str, comp_id: str):
    """Per metal atom, the first-shell donor bond lengths (skips alt confs)."""
    atoms = [(el, xyz) for el, xyz, alt in _atoms(pdb_text) if alt in " A"]
    metals = [xyz for el, xyz in atoms if el == comp_id]
    donors = [xyz for el, xyz in atoms if el in DONOR_ELEMENTS]
    for m in metals:
        d = np.linalg.norm(np.array(donors) - m, axis=1) if donors else np.empty(0)
        shell = d[(d > SHELL[0]) & (d <= SHELL[1])]
        if len(shell) >= MIN_CN:
            yield shell


def build_metal(comp_id: str, label: str, n_structures: int, max_res: float) -> dict:
    bonds: list[float] = []
    counts: list[int] = []
    ids = _search(comp_id, n_structures, max_res)
    for pid in ids:
        try:
            text = urllib.request.urlopen(
                f"https://files.rcsb.org/download/{pid}.pdb", timeout=30
            ).read().decode()
        except Exception:
            continue
        for shell in _sites(text, comp_id):
            bonds.extend(shell.tolist())
            counts.append(len(shell))
    return {
        "metal": label,
        "coordination_number": Counter(counts).most_common(1)[0][0],
        "cn_range": [int(np.percentile(counts, 10)), int(np.percentile(counts, 90))],
        "bond_length_mean": round(statistics.fmean(bonds), 3),
        "bond_length_std": round(statistics.pstdev(bonds), 3),
        "source": f"RCSB PDB ≤{max_res}Å, {len(ids)} structures, "
                  f"{len(counts)} sites / {len(bonds)} bonds, pulled {date.today()}",
    }


def main(n_structures: int = 40, max_res: float = 1.8) -> None:
    table = {}
    for comp_id, label in METALS.items():
        print(f"building {label} from PDB comp {comp_id} ...", flush=True)
        table[label] = build_metal(comp_id, label, n_structures, max_res)
        print(f"  {table[label]}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
