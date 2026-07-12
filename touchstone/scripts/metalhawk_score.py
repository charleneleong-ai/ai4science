"""Score designs with MetalHawk — the open ANN geometry-distortion oracle the
MetalHawkVerifier consumes. Runs on a box with the MetalHawk repo (github.com/vrettasm/
MetalHawk) + its venv (pymol-open-source for site extraction, scikit-learn for the model)
plus gemmi (to normalize inputs — BoltzGen CIFs otherwise crash pymol's metal parsing).

MetalHawk's real pipeline (verified against the released repo):
  1. extract_metal_sites.py --inputdir <dir> --outputdir <dir>  carves a 10 Å sphere per metal
  2. metalhawk.py -f <spheres>/*.pdb [--csd|--no-csd] -o <dir>   classifies each site's geometry
     CLASS, writing a CSV with columns (file, prediction, entropy) — the *geometry class*
     (LIN/TRI/TET/SPL/SQP/TBP/OCT), NOT a CN or probability. CN is implied by the class, and
     confidence is 1 − entropy/ln(7) (ln 7 = max entropy over the 7 classes).

This wraps both and writes {structure_path: {coordination_number, geometry, confidence}} to a
JSON that touchstone.geometry.metalhawk.load_predictions + score_provider consume — the same
file-handoff seam as chai_crosscheck / allmetal3d_crosscheck. The default model is the
MetalPDB-trained one (--no-csd) since designs are proteins. Note: MetalHawk is a *learned*
classifier and goes out-of-distribution on de-novo designs (it over-calls OCT); the
MetalHawkVerifier defers on confident-yet-contradictory calls rather than trusting them.

    conda run -n metalhawk python scripts/metalhawk_score.py \
        'refold_cif/*.cif' --repo ~/MetalHawk --out metalhawk_scores.json
"""

from __future__ import annotations

import csv
import glob
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import gemmi
import typer

# MetalHawk's 7 geometry classes → (coordination number, human-readable geometry).
CLASS = {
    "LIN": (2, "linear"), "TRI": (3, "trigonal planar"), "TET": (4, "tetrahedral"),
    "SPL": (4, "square planar"), "SQP": (5, "square pyramidal"),
    "TBP": (5, "trigonal bipyramidal"), "OCT": (6, "octahedral"),
}
MAX_ENTROPY = math.log(len(CLASS))  # max Shannon entropy over the geometry classes → confidence floor


def normalize(src: str, pdb_dir: Path) -> str:
    """gemmi-normalize a structure (CIF/PDB) to a clean PDB pymol can parse; returns its stem."""
    stem = Path(src).stem
    gemmi.read_structure(src).write_pdb(str(pdb_dir / f"{stem}.pdb"))
    return stem


def extract(repo: Path, pdb_dir: Path, sph_dir: Path) -> None:
    py = sys.executable  # the metalhawk env's interpreter (conda run selects it)
    subprocess.run([py, str(repo / "extract_metal_sites.py"),
                    "--inputdir", str(pdb_dir), "--outputdir", str(sph_dir)], check=True)


def predict(repo: Path, sph_dir: Path, pred_dir: Path, csd: bool) -> None:
    py = sys.executable
    model_flag = "--csd" if csd else "--no-csd"
    subprocess.run([py, str(repo / "metalhawk.py"), "-f", *glob.glob(f"{sph_dir}/*.pdb"),
                    model_flag, "-o", str(pred_dir)], check=True)


def parse(pred_dir: Path) -> dict[str, dict]:
    """MetalHawk CSVs → {design_stem: {coordination_number, geometry, confidence}}. Each sphere
    file is '<stem>_<siteidx>'; per design we keep the most confident (lowest-entropy) site."""
    best: dict[str, tuple[float, dict]] = {}
    for row in (r for f in pred_dir.glob("*.csv") for r in csv.DictReader(f.open())):
        cls = row["prediction"]
        if cls not in CLASS:  # unknown class label ⇒ can't map to a CN, skip
            continue
        stem = row["file"].rsplit("_", 1)[0]
        entropy = float(row["entropy"])
        cn, geometry = CLASS[cls]
        pred = {"coordination_number": cn, "geometry": geometry,
                "confidence": round(max(0.0, 1.0 - entropy / MAX_ENTROPY), 4)}
        if stem not in best or entropy < best[stem][0]:
            best[stem] = (entropy, pred)
    return {stem: pred for stem, (_, pred) in best.items()}


def main(pattern: str, repo: Path = typer.Option(...),
         out: Path = typer.Option(Path("metalhawk_scores.json")),
         csd: bool = typer.Option(False, help="Use the CSD-trained model instead of MetalPDB.")) -> None:
    designs = sorted(glob.glob(pattern))
    work = Path(tempfile.mkdtemp())
    pdb_dir, sph_dir, pred_dir = work / "pdb", work / "spheres", work / "pred"
    pdb_dir.mkdir()
    stem_to_path: dict[str, str] = {}
    for s in designs:
        try:
            stem_to_path[normalize(s, pdb_dir)] = s
        except Exception as e:  # unreadable structure ⇒ no prediction for this design
            print(f"{Path(s).name}: normalize failed ({type(e).__name__}) — skipped", flush=True)
    extract(repo, pdb_dir, sph_dir)
    predict(repo, sph_dir, pred_dir, csd)
    by_stem = parse(pred_dir)

    scores = {stem_to_path[stem]: pred for stem, pred in by_stem.items() if stem in stem_to_path}
    for path, pred in scores.items():
        print(f"{Path(path).name}: {pred['geometry']} CN{pred['coordination_number']} "
              f"(conf {pred['confidence']:.0%})", flush=True)
    out.write_text(json.dumps(scores, indent=2) + "\n")
    print(f"wrote {out}: {len(scores)}/{len(designs)} predictions ({'CSD' if csd else 'MetalPDB'} model)")


if __name__ == "__main__":
    typer.run(main)
