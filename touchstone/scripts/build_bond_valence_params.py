"""Recalibrate the bond-valence R₀ parameters — and tolerances — to the metalloprotein domain.

The shipped R₀ (Brese & O'Keeffe 1991) come from **small-molecule crystals**, where metal–donor
bonds are short and well-ordered. Applied to *protein* sites — longer, lower-CN, entatic — they
under-count systematically: the median real Ni²⁺ metalloprotein site sums to BVS ≈ 1.03 against a
formal valence of 2, so the tier defers ~64% of real sites. Same wrong-domain error as the geometry
prior (docs/experiments/2026-07-13-geometry-prior-wrong-domain.md), in a second always-on tier.

The fix, per metal, from real MetalPDB sites (protein donors, solvent excluded — same domain as the
thing the tier judges):

  1. **R₀**: a single additive shift to the metal's element R₀ values, chosen so the *median*
     real-site BVS lands on the formal valence. Median (not least-squares): a per-element fit
     overfits the wide scatter and produces nonsense on sparse donors (Ni–S: 74 bonds → a negative
     shift). The uniform shift preserves the literature's physical N/O/S offsets and only corrects
     the domain-wide bias.
  2. **Tolerances**: from the real-site |BVS − valence| distribution, so the *data* owns them (as
     cn_range does for geometry), not a hardcoded small-molecule constant. Protein BVS scatter is
     wide — wider for Ni/Co than Cu — so tolerances are per-metal.

The honest consequence, recorded in the source note and the docstring: even recalibrated, BVS is a
**weak** filter in the protein domain (the sum conflates coordination number with bond length across
a heterogeneous domain). Its job is catching egregious over/under-bonding, not fine discrimination —
the ensemble discriminates, not BVS alone.

    uv run python scripts/build_bond_valence_params.py

Only metals with ≥ MIN_SITES real sites are recalibrated; the rest keep their literature R₀ with a
flag, so a bundled parameter is never a fit to a handful of sites.
"""

from __future__ import annotations

import json
import math
import statistics
import urllib.parse
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import typer

from touchstone.geometry.parse import DONOR_ELEMENTS, SOLVENT_RESIDUES

OUT = Path(__file__).parent.parent / "src" / "touchstone" / "data" / "bond_valence_params.json"
API = "http://metalpdb.cerm.unifi.it/api?query="
SHELL = (1.0, 2.8)
MIN_CN = 3
MIN_SITES = 30  # too few real sites to trust a recalibration — keep the literature value, flagged

# Literature R₀ (Brese & O'Keeffe 1991), the small-molecule starting point, plus the metal → element
# oxidation state and the MetalPDB element symbol. b is the universal bond-valence constant.
B = 0.37
LITERATURE = {
    "Ni2+": {"symbol": "Ni", "valence": 2, "r0": {"N": 1.647, "O": 1.654, "S": 2.04}},
    "Cu2+": {"symbol": "Cu", "valence": 2, "r0": {"N": 1.713, "O": 1.679, "S": 2.054}},
    "Co2+": {"symbol": "Co", "valence": 2, "r0": {"N": 1.658, "O": 1.692, "S": 2.03}},
}


def fetch(metal_symbol: str) -> list[dict]:
    q = urllib.parse.quote(f"metal:{metal_symbol},representative:TRUE")
    with urllib.request.urlopen(API + q, timeout=120) as resp:
        data = json.load(resp)
    return data if isinstance(data, list) else data.get("sites", data.get("results", []))


def donor_bonds(records: list[dict], metal_symbol: str) -> list[list[tuple[str, float]]]:
    """Per real metal site, its (element, distance) first-shell protein donors."""
    sites = []
    for rec in records:
        for m in rec.get("metals", []):
            if (m.get("symbol") or "").capitalize() != metal_symbol:
                continue
            bonds = [
                ((d.get("symbol") or "").upper(), d["distance"])
                for lig in m.get("ligands", [])
                if (lig.get("residue") or "").upper() not in SOLVENT_RESIDUES
                for d in lig.get("donors", [])
                if (d.get("symbol") or "").upper() in DONOR_ELEMENTS
                and isinstance(d.get("distance"), (int, float))
                and SHELL[0] < d["distance"] <= SHELL[1]
            ]
            if len(bonds) >= MIN_CN:
                sites.append(bonds)
    return sites


def bvs(bonds: list[tuple[str, float]], r0: dict[str, float]) -> float:
    return sum(math.exp((r0[e] - d) / B) for e, d in bonds)


def recalibrate(sites: list[list[tuple[str, float]]], r0: dict[str, float], valence: int) -> dict:
    """Median-centred R₀ shift + data-driven tolerances for one metal."""
    median_bvs = statistics.median(bvs(s, r0) for s in sites)
    shift = B * math.log(valence / median_bvs)  # exp(shift/b) rescales the median onto the valence
    r0_new = {e: round(v + shift, 3) for e, v in r0.items()}
    devs = [abs(bvs(s, r0_new) - valence) for s in sites]
    return {
        "N": r0_new["N"], "O": r0_new["O"], "S": r0_new["S"],
        # the data owns the tolerance: ~68% of real sites within trust, ~90% not deferred
        "trust_tol": round(float(np.percentile(devs, 68)), 2),
        "ood_tol": round(float(np.percentile(devs, 90)), 2),
    }


def main(max_sites: int = 3000) -> None:
    metals = {}
    for label, spec in LITERATURE.items():
        print(f"recalibrating {label} from MetalPDB element {spec['symbol']} ...", flush=True)
        try:
            sites = donor_bonds(fetch(spec["symbol"]), spec["symbol"])[:max_sites]
        except Exception as e:
            print(f"  {label}: fetch failed ({type(e).__name__}) — keeping literature R0", flush=True)
            metals[label] = {**spec["r0"], "domain": "small-molecule (Brese & O'Keeffe 1991)"}
            continue
        if len(sites) < MIN_SITES:
            print(f"  {label}: only {len(sites)} sites (< {MIN_SITES}) — keeping literature R0", flush=True)
            metals[label] = {**spec["r0"], "domain": "small-molecule (Brese & O'Keeffe 1991)"}
            continue
        entry = recalibrate(sites, spec["r0"], spec["valence"])
        entry["domain"] = f"metalloprotein (MetalPDB, {len(sites)} sites, protein donors)"
        metals[label] = entry
        old_defer = sum(abs(bvs(s, spec["r0"]) - spec["valence"]) > 0.8 for s in sites) / len(sites)
        new_defer = sum(abs(bvs(s, entry) - spec["valence"]) > entry["ood_tol"] for s in sites) / len(sites)
        print(f"  {entry}   real-site defer {old_defer:.0%} → {new_defer:.0%}", flush=True)

    if not any(m.get("domain", "").startswith("metalloprotein") for m in metals.values()):
        raise SystemExit("no metal recalibrated (offline? API down?) — refusing to overwrite params")

    out = {
        "b": B,
        "metals": metals,
        "source": (
            f"R0 recalibrated to the metalloprotein domain from MetalPDB, pulled {date.today()} "
            "(median-centred shift on Brese & O'Keeffe 1991 R0; per-metal tolerances from the "
            "real-site BVS scatter). BVS is a weak filter in this domain — it catches egregious "
            "over/under-bonding, not fine metal discrimination. See "
            "docs/experiments/2026-07-17-bond-valence-wrong-domain.md"
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    typer.run(main)
