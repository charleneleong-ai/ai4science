"""Score designs with MetalHawk — the open ANN geometry-distortion oracle the
MetalHawkVerifier consumes. Runs on a box with the MetalHawk repo (github.com/vrettasm/
MetalHawk) + its venv (pymol-open-source for the site extraction, torch for the ANN).

MetalHawk's pipeline: extract_metal_sites.py carves a 10 Å sphere around each metal, then
metalhawk.py classifies the site's coordination number + geometry with a confidence. This
wraps both and writes {structure_path: {coordination_number, geometry, confidence}} to a
JSON that touchstone.geometry.metalhawk.score_provider loads — the same file-handoff seam
as chai_crosscheck / allmetal3d_crosscheck.

    conda run -n metalhawk python scripts/metalhawk_score.py \
        'bg_motif_pdbs/*.pdb' --repo ~/MetalHawk --out metalhawk_scores.json

The output-parsing (_parse) targets MetalHawk's released CSV columns (file, CN, geometry,
probability); adjust it to your installed version if the columns differ.
"""

from __future__ import annotations

import csv
import glob
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import typer


def _run_metalhawk(structure: Path, repo: Path, work: Path) -> Path:
    """extract_metal_sites → metalhawk; returns the dir holding MetalHawk's output."""
    spheres, out = work / "spheres", work / "pred"
    spheres.mkdir(parents=True, exist_ok=True)
    py = sys.executable  # the metalhawk env's interpreter (conda run selects it)
    subprocess.run([py, str(repo / "extract_metal_sites.py"), "-i", str(structure), "-o", str(spheres)], check=True)
    subprocess.run([py, str(repo / "metalhawk.py"), "-f", *glob.glob(f"{spheres}/*.pdb"), "-o", str(out)], check=True)
    return out


def _parse(out: Path) -> tuple[int, str, float] | None:
    """(coordination_number, geometry, confidence) from MetalHawk's top prediction."""
    rows = [r for f in out.glob("*.csv") for r in csv.DictReader(f.open())]
    if not rows:
        return None
    top = max(rows, key=lambda r: float(r.get("probability", r.get("confidence", 0)) or 0))
    cn = int(float(top.get("CN", top.get("coordination_number", 0))))
    geometry = top.get("geometry", top.get("class", "unknown"))
    confidence = float(top.get("probability", top.get("confidence", 0)) or 0)
    return cn, geometry, confidence


def main(pattern: str, repo: Path = typer.Option(...), out: Path = typer.Option(Path("metalhawk_scores.json"))) -> None:
    scores: dict[str, dict] = {}
    for s in sorted(glob.glob(pattern)):
        work = Path(tempfile.mkdtemp())
        try:
            pred = _parse(_run_metalhawk(Path(s), repo, work))
        except Exception as e:  # extraction / inference failure ⇒ no prediction for this design
            print(f"{Path(s).name}: MetalHawk failed ({type(e).__name__}) — skipped", flush=True)
            continue
        if pred is None:
            print(f"{Path(s).name}: no MetalHawk prediction — skipped", flush=True)
            continue
        cn, geometry, confidence = pred
        scores[s] = {"coordination_number": cn, "geometry": geometry, "confidence": confidence}
        print(f"{Path(s).name}: {geometry} CN{cn} (conf {confidence:.0%})", flush=True)
    out.write_text(json.dumps(scores, indent=2) + "\n")
    print(f"wrote {out}: {len(scores)} predictions")


if __name__ == "__main__":
    typer.run(main)
