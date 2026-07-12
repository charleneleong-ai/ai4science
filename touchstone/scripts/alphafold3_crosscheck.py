"""Co-fold a design's sequence + metal with AlphaFold3 — a fourth, independent
prediction for the CofoldCrossCheck verifier. Runs on a GPU box with AF3's inference
code (google-deepmind/alphafold3) + its (gated) model weights.

AF3 is a distinct model from Boltz-2 / Chai-1, so it is a genuinely independent re-fold:
does the generator's predicted metal site reproduce when AF3 folds the sequence from
scratch? It places the metal as a CCD ligand entry (NI / CU / CO / ZN / FE / MN / …). The
output mmCIF is parsed by touchstone's `cif_provider` / `coordination_site` identically to
any other structure, so AF3 drops in as a provider with no verifier change.

Weights are gated (request access from DeepMind): point `--model-dir` at them and `--db-dir`
at the MSA databases. The AF3 input-JSON `version` may need bumping to match your AF3 release.

    conda run -n af3 python scripts/alphafold3_crosscheck.py \
        --seq MKLVINGKTLKGEITV... --name design_4 --out af3_out/ --metal NI \
        --model-dir $AF3_MODELS --db-dir $AF3_DBS
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import typer

CCD = {"Ni2+": "NI", "Cu2+": "CU", "Co2+": "CO", "Zn2+": "ZN", "Fe2+": "FE", "Mn2+": "MN"}


def main(
    seq: str = typer.Option(..., help="protein sequence (one-letter)"),
    out: Path = typer.Option(..., help="output directory for the predicted CIF"),
    name: str = "design",
    metal: str = typer.Option("NI", help="metal ion as a PDB CCD code (NI / CU / CO / ZN / …) or a touchstone label (Ni2+)"),
    model_dir: Path = typer.Option(..., help="AF3 model-weights dir (gated — request from DeepMind)"),
    db_dir: Path = typer.Option(..., help="AF3 sequence-database dir (for the MSA pipeline)"),
    run_alphafold: Path = typer.Option("run_alphafold.py", help="path to AF3's run_alphafold.py"),
    seed: int = 42,
) -> None:
    ccd = CCD.get(metal, metal).upper()  # accept either a CCD code or a touchstone metal label
    out.mkdir(parents=True, exist_ok=True)
    spec = {
        "name": name,
        "modelSeeds": [seed],
        "sequences": [
            {"protein": {"id": "A", "sequence": seq}},
            {"ligand": {"id": "B", "ccdCodes": [ccd]}},
        ],
        "dialect": "alphafold3",
        "version": 2,
    }
    json_path = out / f"{name}.json"
    json_path.write_text(json.dumps(spec, indent=2))

    subprocess.run(
        ["python", str(run_alphafold), f"--json_path={json_path}",
         f"--model_dir={model_dir}", f"--db_dir={db_dir}", f"--output_dir={out}"],
        check=True,
    )
    cifs = sorted(out.glob("**/*model*.cif")) or sorted(out.glob("**/*.cif"))
    print(f"wrote {len(cifs)} prediction(s) to {out}  (metal ligand: {ccd})")
    for c in cifs:
        print(f"  {c}")


if __name__ == "__main__":
    typer.run(main)
