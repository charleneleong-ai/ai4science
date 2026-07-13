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
from collections.abc import Iterator
from datetime import date
from pathlib import Path

import numpy as np
import typer

from touchstone.geometry.parse import DONOR_ELEMENTS, SOLVENT_RESIDUES

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
# Don't ship a prior we can't support. The precious metals are sparse in the PDB (Au³⁺: 4 sites,
# Pd²⁺: 6) — a mean±std from a handful of sites is noise with a decimal point, and a bundled
# reference is trusted silently by every verdict. Below this floor, emit nothing and let the
# metal report "no reference geometry" rather than a fabricated one.
MIN_SITES = 30
# Solvent donors are excluded: a *designed* structure carries no waters, so touchstone can only ever
# measure protein donors. Counting MetalPDB's waters would compare a design against donors it
# structurally cannot present — Ni/Co/Mn sit octahedrally with ~2 aqua ligands, so a water-inclusive
# prior reports modal CN 6 and marks a perfectly good 4-protein-donor site (open face, water-fillable
# in solution) as two donors short. The prior must live in the same domain as the thing it judges,
# so the residue set is imported from the parser rather than duplicated here.


def fetch(metal_symbol: str) -> list[dict]:
    """MetalPDB representative sites for a metal element."""
    q = urllib.parse.quote(f"metal:{metal_symbol},representative:TRUE")
    with urllib.request.urlopen(API + q, timeout=120) as resp:
        data = json.load(resp)
    return data if isinstance(data, list) else data.get("sites", data.get("results", []))


def donor_shells(records: list[dict], metal_symbol: str) -> Iterator[list[float]]:
    """Per metal site, the first-shell N/O/S donor distances (Å), protein donors only (see SOLVENT_RESIDUES).

    MIN_CN is a domain filter, not just noise rejection: 58% of MetalPDB Ni sites hold the metal with
    only 1-2 protein contacts and a shell of water. Those are surface/adventitious ions, not the
    buried sites a design targets. Dropping them moves the bond mean by 0.013 Å (2.203 → 2.190), so
    the selection costs almost nothing — but it does mean `cn_range`'s lower bound is this floor
    rather than a free observation.
    """
    for rec in records:
        for m in rec.get("metals", []):
            if (m.get("symbol") or "").capitalize() != metal_symbol:
                continue
            dists = [
                d["distance"]
                for lig in m.get("ligands", [])
                if (lig.get("residue") or "").upper() not in SOLVENT_RESIDUES
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
        # floor/ceil, never int()/round(): both truncate *inward*, narrowing the trust gate. A true
        # 90th-pct CN of 5.9 must ship a ceiling of 6 (int() gives 5, rejecting real 6-coordinate
        # sites); a 10th-pct of 3.6 must ship a floor of 3 (round() gives 4, rejecting real CN-3).
        "cn_range": [int(np.floor(np.percentile(counts, 10))), int(np.ceil(np.percentile(counts, 90)))],
        "bond_length_mean": round(statistics.fmean(bonds), 3),
        "bond_length_std": round(statistics.pstdev(bonds), 3),
        "source": source,
    }


def main(max_sites: int = 2000) -> None:
    table = {}
    for symbol, label in METALS.items():
        print(f"building {label} from MetalPDB element {symbol} ...", flush=True)
        try:
            shells = list(donor_shells(fetch(symbol), symbol))[:max_sites]
        except Exception as e:  # a per-metal network/parse failure shouldn't sink the batch
            print(f"  {label}: fetch failed ({type(e).__name__}) — skipped", flush=True)
            continue
        if len(shells) < MIN_SITES:  # too sparse to support a prior — ship nothing, not noise
            print(f"  {label}: only {len(shells)} sites (< {MIN_SITES}) — skipped", flush=True)
            continue
        note = f"MetalPDB representative sites (protein donors, solvent excluded), {len(shells)} sites, pulled {date.today()}"
        table[label] = aggregate(shells, label, note)
        print(f"  {table[label]}", flush=True)
    if not table:
        # never leave an empty table on disk: best_reference() selects on file *existence*, so an
        # offline run would install a MetalPDB prior that KeyErrors on every metal — the geometry
        # tier dies silently. Absent beats fabricated; present-but-empty is worse than either.
        raise SystemExit("no metal met the site floor (offline? API down?) — refusing to overwrite the prior")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
