"""Build touchstone's open coordination-motif precedent table from MetalPDB.

The precedent tier (touchstone.geometry.precedent) asks: has this metal–donor motif been seen
in nature? CSD-CrossMiner answers it against the CSD (licence-gated); this builds the open,
licence-free analog from MetalPDB (metalpdb.cerm.unifi.it) — the same source as the geometry
reference (build_metalpdb_reference.py). For each metal it pulls representative PDB sites, labels
each by its first-shell donor motif (e.g. 'Ni-N3O2'), and writes {motif: n_sites} to
src/touchstone/data/metalpdb_precedents.json, which precedent.metalpdb_precedent_search reads.

    uv run python scripts/build_metalpdb_precedents.py
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

import typer

from touchstone.geometry.parse import DONOR_ELEMENTS, SOLVENT_RESIDUES

OUT = Path(__file__).parent.parent / "src" / "touchstone" / "data" / "metalpdb_precedents.json"
# element → verifier label; the metal element symbol is also the motif prefix (matches
# precedent.motif, which keys on element_symbol(metal)).
METALS = {
    "Ni": "Ni2+", "Cu": "Cu2+", "Co": "Co2+",  # metal-recovery targets
    "Zn": "Zn2+", "Fe": "Fe3+", "Mn": "Mn2+",  # common biological transition metals
    "Ca": "Ca2+", "Mg": "Mg2+",                # alkaline-earth / hard cations
    "Pd": "Pd2+", "Pt": "Pt2+", "Au": "Au3+",  # precious (e-waste recovery)
}
API = "http://metalpdb.cerm.unifi.it/api?query="
SHELL = (1.0, 2.8)  # first-shell distance window, Å
MIN_CN = 3  # ignore adventitious ions with < 3 contacts


def fetch(metal_symbol: str) -> list[dict]:
    """MetalPDB representative sites for a metal element."""
    q = urllib.parse.quote(f"metal:{metal_symbol},representative:TRUE")
    with urllib.request.urlopen(API + q, timeout=120) as resp:
        data = json.load(resp)
    return data if isinstance(data, list) else data.get("sites", data.get("results", []))


def site_motifs(records: list[dict], metal_symbol: str) -> Iterator[str]:
    """Per metal site, its first-shell donor-composition motif (e.g. 'Ni-N3O2'), protein donors only.

    Solvent is excluded for the same reason as in the geometry prior: a design carries no waters, so
    a motif it can never present must not compete for precedent counts. Counting them buries the real
    protein motifs under aqua-padded variants — with waters in, `Ni-O6` (hexaaqua nickel, an ion
    sitting in solvent) scores 22 hits while `Ni-S4`, the NiFe-hydrogenase Cys4 site, scores 3 and
    fails the min_hits gate.
    """
    for rec in records:
        for m in rec.get("metals", []):
            if (m.get("symbol") or "").capitalize() != metal_symbol:
                continue
            donors = [
                (d.get("symbol") or "").upper()
                for lig in m.get("ligands", [])
                if (lig.get("residue") or "").upper() not in SOLVENT_RESIDUES
                for d in lig.get("donors", [])
                if (d.get("symbol") or "").upper() in DONOR_ELEMENTS
                and isinstance(d.get("distance"), (int, float))
                and SHELL[0] < d["distance"] <= SHELL[1]
            ]
            if len(donors) >= MIN_CN:
                comp = "".join(f"{el}{n}" for el, n in sorted(Counter(donors).items()))
                yield f"{metal_symbol}-{comp}"


def main(max_sites: int = 5000) -> None:
    table: Counter[str] = Counter()
    for symbol, label in METALS.items():
        print(f"mining {label} ({symbol}) precedents from MetalPDB ...", flush=True)
        try:
            motifs = list(site_motifs(fetch(symbol), symbol))[:max_sites]
        except Exception as e:  # a per-metal network/parse failure shouldn't sink the batch
            print(f"  {label}: fetch failed ({type(e).__name__}) — skipped", flush=True)
            continue
        table.update(motifs)
        print(f"  {label}: {len(motifs)} sites, {len(set(motifs))} distinct motifs", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(table.most_common()), indent=2) + "\n")
    print(f"wrote {OUT}: {len(table)} motifs across {sum(table.values())} sites, pulled {date.today()}")


if __name__ == "__main__":
    typer.run(main)
