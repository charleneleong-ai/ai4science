"""Build touchstone's reference metal geometry from MetalPDB — the open, licence-free
metalloprotein analog of the CSD prior.

MetalPDB (metalpdb.cerm.unifi.it) curates every metal site in the PDB with its
coordinating donors, distances, and coordination number. It is metalloprotein-specific
(unlike the CSD's small-molecule crystals) and needs no licence, so it is the natural
default reference for designed metal-binding *proteins*. Pulls each metal's representative
sites, aggregates first-shell donor distances + modal CN, and writes the empirical
mean/std + CN range to src/touchstone/data/metalpdb_reference.json — schema identical to
pdb/csd_reference.json, so MetalPDBReference drops in behind the same geometry() interface.

    uv run python scripts/build_metalpdb_reference.py

REST API (see protti::fetch_metal_pdb): GET metalpdb.cerm.unifi.it/api?query=metal:Ni,...
returns sites nested metals → ligands → donors → {symbol, distance}.
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

OUT = Path(__file__).parent.parent / "src" / "touchstone" / "data" / "metalpdb_reference.json"
# MetalPDB element symbol → verifier label. Open/licence-free, so the panel spans the
# recovery targets and the common biological + precious metals (element-level pull; one
# representative oxidation state per element — Fe mixes Fe²⁺/Fe³⁺ sites).
METALS = {
    "Ni": "Ni2+", "Cu": "Cu2+", "Co": "Co2+",  # metal-recovery targets
    "Zn": "Zn2+", "Fe": "Fe3+", "Mn": "Mn2+",  # common biological transition metals
    "Ca": "Ca2+", "Mg": "Mg2+",                # alkaline-earth / hard cations
    "Pd": "Pd2+", "Pt": "Pt2+", "Au": "Au3+",  # precious (e-waste recovery)
}
API = "http://metalpdb.cerm.unifi.it/api?query="
SHELL = (1.0, 2.8)  # first-shell distance window, Å
MIN_CN = 3  # ignore adventitious ions with < 3 contacts


def _fetch(metal_symbol: str) -> list[dict]:
    """MetalPDB representative sites for a metal element."""
    q = urllib.parse.quote(f"metal:{metal_symbol},representative:TRUE")
    with urllib.request.urlopen(API + q, timeout=120) as resp:
        data = json.load(resp)
    return data if isinstance(data, list) else data.get("sites", data.get("results", []))


def _shells(records: list[dict], metal_symbol: str):
    """Per metal site, the first-shell N/O/S donor distances (Å)."""
    for rec in records:
        for m in rec.get("metals", []):
            if (m.get("symbol") or "").capitalize() != metal_symbol:
                continue
            dists = [
                d["distance"]
                for lig in m.get("ligands", [])
                for d in lig.get("donors", [])
                if (d.get("symbol") or "").upper() in DONOR_ELEMENTS and isinstance(d.get("distance"), (int, float))
            ]
            shell = [x for x in dists if SHELL[0] < x <= SHELL[1]]
            if len(shell) >= MIN_CN:
                yield shell


def aggregate(shells: list[list[float]], label: str, source: str) -> dict:
    """Empirical geometry (mean/std bond, modal CN, CN range) from per-site donor shells —
    pure, so the JSON schema is testable without hitting the network."""
    bonds = [b for shell in shells for b in shell]
    counts = [len(shell) for shell in shells]
    return {
        "metal": label,
        "coordination_number": Counter(counts).most_common(1)[0][0],
        "cn_range": [int(np.percentile(counts, 10)), int(np.percentile(counts, 90))],
        "bond_length_mean": round(statistics.fmean(bonds), 3),
        "bond_length_std": round(statistics.pstdev(bonds), 3),
        "source": source,
    }


def main(max_sites: int = 2000) -> None:
    table = {}
    for symbol, label in METALS.items():
        print(f"building {label} from MetalPDB element {symbol} ...", flush=True)
        try:
            shells = list(_shells(_fetch(symbol), symbol))[:max_sites]
        except Exception as e:  # a per-metal network/parse failure shouldn't sink the batch
            print(f"  {label}: fetch failed ({type(e).__name__}) — skipped", flush=True)
            continue
        if not shells:  # sparse metals (some precious ones) may have too few PDB sites
            print(f"  {label}: no qualifying MetalPDB sites — skipped", flush=True)
            continue
        note = f"MetalPDB representative sites, {len(shells)} sites, pulled {date.today()}"
        table[label] = aggregate(shells, label, note)
        print(f"  {table[label]}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
