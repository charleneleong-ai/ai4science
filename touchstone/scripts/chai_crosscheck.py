"""Co-fold a design's sequence + metal with Chai-1 — the independent prediction
the CofoldCrossCheck verifier scores against. Runs on a GPU box with chai_lab.

Chai-1 is an open (Apache-2.0) AlphaFold3-class co-folder, so it is genuinely
independent of Boltz-2: a second predictor, not a re-run of the first.

    conda run -n chai python scripts/chai_crosscheck.py \
        --seq MKLVINGKTLKGEITV... --name design_4 --out chai_out/ --metal "[Ni+2]"
"""

from __future__ import annotations

from pathlib import Path

import typer


def main(
    seq: str = typer.Option(..., help="protein sequence (one-letter)"),
    out: Path = typer.Option(..., help="output directory for the predicted CIF"),
    name: str = "design",
    metal: str = typer.Option("[Ni+2]", help="metal ion as SMILES (Chai ligand entry)"),
    seed: int = 42,
) -> None:
    from chai_lab.chai1 import run_inference

    out.mkdir(parents=True, exist_ok=True)
    fasta = out / f"{name}.fasta"
    fasta.write_text(f">protein|name={name}\n{seq}\n>ligand|name=METAL\n{metal}\n")
    pred_dir = out / f"{name}_pred"  # Chai requires a fresh/empty output dir

    candidates = run_inference(
        fasta_file=fasta,
        output_dir=pred_dir,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=seed,
        device="cuda",
        use_esm_embeddings=True,
    )
    cifs = sorted(pred_dir.glob("*.cif"))
    print(f"wrote {len(cifs)} prediction(s) to {out}")
    for c in cifs:
        print(f"  {c}")
    if candidates is not None and getattr(candidates, "ranking_data", None):
        print("ranked by Chai aggregate score (best first)")


if __name__ == "__main__":
    typer.run(main)
