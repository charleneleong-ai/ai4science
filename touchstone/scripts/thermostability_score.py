"""Predict a design's melting temperature for the ThermostabilityVerifier.

Real predictor: **TemStaPro** (Pudžiuvelytė et al. 2024) — ProtT5 embeddings + an
ensemble of binary classifiers for thermostability at 40/45/50/55/60/65 °C. Its
`--mean-output` TSV gives, per sequence, a "stable above T?" call at each threshold;
`tm_from_thresholds` collapses that monotonic profile to one Tm estimate (the highest
threshold the protein is predicted to stay folded above) — the single number the
`ThermostabilityVerifier` / `tm_provider` consume. So TemStaPro's native classification
drives the verdict directly, just summarised to a Tm.

Runs on a GPU box in TemStaPro's own conda env (ProtT5 download is automatic to `-d`):

    git clone https://github.com/ievapudz/TemStaPro.git    # (run via `!` — deny-listed for the agent)
    conda env create -f TemStaPro/environment_GPU.yml       # env temstapro_env_GPU
    conda run -n temstapro_env_GPU python scripts/thermostability_score.py \
        --fasta designs.fasta --temstapro-dir TemStaPro --out tm.json

`tm.json` is a {sequence-id: Tm °C} map; pair it with the sequences and feed
`tm_provider({sequence: Tm})`. `--placeholder` runs the no-model heuristic for a smoke test.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import typer

THRESHOLDS = (40, 45, 50, 55, 60, 65)  # TemStaPro's binary-classifier temperatures, °C


def tm_from_thresholds(passed: dict[int, bool]) -> float:
    """Collapse TemStaPro's per-threshold 'stable above T?' calls into one Tm estimate:
    the highest threshold the protein is predicted to stay folded above. Below the lowest
    threshold ⇒ ~37 °C (mesophilic/unstable); above the highest ⇒ ~70 °C (thermophilic —
    the verifier's interest saturates past its trust cutoff anyway)."""
    highest = max((t for t, ok in passed.items() if ok), default=None)
    if highest is None:
        return 37.0
    if highest == THRESHOLDS[-1]:
        return 70.0
    return float(highest)


def placeholder_tm(seq: str) -> float:
    """No-model smoke-test Tm (°C): charged/polar fraction proxy. NOT for quantitative use."""
    charged = sum(a in "DEKR" for a in seq) / max(len(seq), 1)
    polar = sum(a in "STNQHY" for a in seq) / max(len(seq), 1)
    return min(max(35.0 + 60.0 * charged + 20.0 * polar, 0.0), 110.0)


def binary_columns(fieldnames: list[str]) -> dict[int, str]:
    """Map each threshold to its binary-prediction column in the mean-output TSV. TemStaPro
    names the binary calls by the threshold number; prefer an exact `str(T)` match, else the
    column that mentions T and isn't the 'raw'/probability one. Validate against a real run."""
    cols: dict[int, str] = {}
    for t in THRESHOLDS:
        exact = [c for c in fieldnames if c.strip() == str(t)]
        fuzzy = [c for c in fieldnames if str(t) in c and "raw" not in c.lower()]
        if exact or fuzzy:
            cols[t] = (exact or fuzzy)[0]
    return cols


def run_temstapro(fasta: Path, temstapro_dir: Path, prottrans_dir: Path, cache_dir: Path) -> dict[str, float]:
    """Run TemStaPro `--mean-output` over a FASTA and return {sequence-id: Tm °C}."""
    out_tsv = cache_dir / "mean_output.tsv"
    cache_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(temstapro_dir / "temstapro"), "-f", str(fasta), "-d", str(prottrans_dir),
         "-e", str(cache_dir), "--mean-output", str(out_tsv)],
        check=True,
    )
    with out_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = binary_columns(reader.fieldnames or [])
        if len(cols) != len(THRESHOLDS):
            raise ValueError(f"could not map all thresholds to TSV columns: got {cols} from {reader.fieldnames}")
        id_col = (reader.fieldnames or ["sequence"])[0]
        return {
            row[id_col]: tm_from_thresholds({t: str(row[c]).strip() in ("1", "1.0", "True") for t, c in cols.items()})
            for row in reader
        }


def main(
    fasta: Path = typer.Option(None, help="FASTA of designs (batch — the efficient path)"),
    temstapro_dir: Path = typer.Option(None, help="cloned TemStaPro repo"),
    prottrans_dir: Path = typer.Option("./ProtTrans", help="ProtT5 cache dir (auto-downloaded)"),
    cache_dir: Path = typer.Option("./temstapro_cache", help="embeddings + output cache"),
    out: Path = typer.Option(None, help="write {sequence-id: Tm} JSON here"),
    seq: str = typer.Option(None, help="single sequence (use with --placeholder for a smoke test)"),
    placeholder: bool = typer.Option(False, help="skip TemStaPro, use the no-model heuristic"),
) -> None:
    if placeholder or (seq and not temstapro_dir):
        if not seq:
            raise typer.BadParameter("--placeholder needs --seq")
        print(f"tm={placeholder_tm(seq):.1f}")
        return
    if not (fasta and temstapro_dir):
        raise typer.BadParameter("real prediction needs --fasta and --temstapro-dir")
    tms = run_temstapro(fasta, temstapro_dir, prottrans_dir, cache_dir)
    if out:
        out.write_text(json.dumps(tms, indent=2) + "\n")
        print(f"wrote {out} ({len(tms)} sequences)")
    else:
        for sid, tm in tms.items():
            print(f"{sid}\ttm={tm:.1f}")


if __name__ == "__main__":
    typer.run(main)
