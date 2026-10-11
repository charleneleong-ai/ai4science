"""Build the full His-switch campaign: every design with switchable positions, parent + variants.

Rationale from the measured data:
  - 200 single substitutions scanned, 18 ideal (pKa 6.3-7.2) and 17 usable (6.0-7.6).
  - No structural descriptor predicts which positions work (burial r=-0.37, carboxylate count
    r=+0.21, cation count the wrong sign), so positions come from the pKa scan, not a rule.
  - On the one design tested so far, ALL FOUR His variants improved binder-target ipTM over the
    parent (0.288 -> 0.366-0.444 mean). That is n=1, so this replicates it across designs.

Emits Boltz co-folding specs for each parent and its variants, so parent and variant are
always measured under identical conditions.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from pathlib import Path
from typing import Any

import typer
from rich.console import Console

from tnf_binder.rosetta import HisMutator
from tnf_binder.specs import cofold_spec
from tnf_binder.tables import read_rows, write_rows
from tnf_binder.target import HUMAN

app = typer.Typer(add_completion=False)
console = Console()

MANIFEST_COLUMNS = ["spec", "design", "kind", "positions", "pkas", "sequence"]
# The middle of the switch window: variants are ordered by how close their pKa sits to it.
WINDOW_CENTRE = 6.7

ScanSite = tuple[str, int, str, float, str]


def switchable_by_design(scans: list[Path]) -> dict[str, list[ScanSite]]:
    """design -> switchable positions, ideal first, then nearest the middle of the window.

    Positions come from the pKa scan, never from a structural rule: burial correlates only
    r=-0.37 with predicted pKa, and carboxylate and cation counts fare worse still.
    """
    found: dict[str, list[ScanSite]] = defaultdict(list)
    for path in scans:
        for row in read_rows(path):
            if row.get("switchable", "").strip().lower() == "true":
                found[row["parent"]].append((row["chain"], int(row["resnum"]), row["wt"],
                                             float(row["pka"]), row["verdict"]))
    for sites in found.values():
        sites.sort(key=lambda site: (site[4] != "ideal", abs(site[3] - WINDOW_CENTRE)))
    return dict(found)


def find_pdb(stem: str, directories: list[Path]) -> Path | None:
    return next((d / f"{stem}.pdb" for d in directories if (d / f"{stem}.pdb").exists()), None)


def variant_combos(sites: list[ScanSite], max_his: int, cap: int) -> list[tuple[ScanSite, ...]]:
    combos: list[tuple[ScanSite, ...]] = []
    for k in range(1, min(max_his, len(sites)) + 1):
        combos += list(itertools.combinations(sites, k))
    return combos[:cap]


def variant_rows(mutator: HisMutator, parent_pose: Any, design: str,
                 combos: list[tuple[ScanSite, ...]], binder_chain: str,
                 out_dir: Path) -> list[dict[str, str]]:
    rows = []
    for combo in combos:
        pose = parent_pose.clone()
        if not mutator.mutate(pose, [(chain, num) for chain, num, *_ in combo]):
            continue
        tag = "_".join(f"{chain}{num}H" for chain, num, *_ in combo)
        name = f"V__{design}__{tag}"
        sequence = mutator.chain_sequence(pose, binder_chain)
        (out_dir / f"{name}.yaml").write_text(cofold_spec(HUMAN, binder=sequence))
        rows.append({"spec": name, "design": design, "kind": "variant",
                     "positions": " ".join(f"{chain}{num}" for chain, num, *_ in combo),
                     "pkas": " ".join(f"{chain}{num}={pka:.2f}"
                                      for chain, num, _wt, pka, _v in combo),
                     "sequence": sequence})
    return rows


@app.command()
def main(
    scans: list[Path] = typer.Argument(..., help="pka scan CSVs from his_pka_scan.py"),
    pdb_dirs: list[Path] = typer.Option(..., "--pdb-dir", help="where the parent PDBs live"),
    out_dir: Path = typer.Option(..., "--out-dir", help="where the Boltz specs go"),
    manifest: Path = typer.Option(..., "--manifest", help="CSV mapping spec -> design + pKas"),
    max_his: int = typer.Option(3, "--max-his"),
    max_variants: int = typer.Option(4, "--max-variants", help="variants per design"),
) -> None:
    """Write a parent spec and up to `max_variants` His-variant specs per scanned design."""
    by_design = switchable_by_design(scans)
    console.print(f"[bold]{len(by_design)}[/] designs carry switchable positions")
    out_dir.mkdir(parents=True, exist_ok=True)

    mutator = HisMutator()
    rows: list[dict[str, str]] = []
    for design, sites in sorted(by_design.items()):
        source = find_pdb(design, pdb_dirs)
        if source is None:
            console.print(f"  {design}: no pdb found, skipping")
            continue
        parent_pose = mutator.load(source)
        binder_chain = mutator.binder_chain(parent_pose)
        parent_sequence = mutator.chain_sequence(parent_pose, binder_chain)
        (out_dir / f"P__{design}.yaml").write_text(cofold_spec(HUMAN, binder=parent_sequence))
        rows.append({"spec": f"P__{design}", "design": design, "kind": "parent",
                     "positions": "", "pkas": "", "sequence": parent_sequence})

        on_binder = [s for s in sites if s[0] == binder_chain][:4]
        combos = variant_combos(on_binder, max_his, max_variants)
        rows += variant_rows(mutator, parent_pose, design, combos, binder_chain, out_dir)
        console.print(f"  {design}: binder chain {binder_chain}, "
                      f"{len(on_binder)} switchable -> {len(combos)} variants")

    write_rows(manifest, rows, MANIFEST_COLUMNS)
    console.print(f"\nwrote [bold]{len(rows)}[/] specs -> [cyan]{out_dir}[/] "
                  f"(manifest [cyan]{manifest}[/])")


if __name__ == "__main__":
    app()
