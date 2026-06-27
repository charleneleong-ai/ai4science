"""Predict a design's melting temperature for the ThermostabilityVerifier. Runs on a
GPU box with a sequence Tm regressor.

This is a SCAFFOLD: wire in a real predictor — TemStaPro or DeepSTABp (sequence → Tm,
ProtTrans-based) or ThermoMPNN (structure-aware ΔΔG) — at `_predict_tm`. The placeholder
returns a crude charge/length heuristic so the pipeline runs end to end; replace before
quantitative use.

    conda run -n thermo python scripts/thermostability_score.py --seq MKLV...
"""

from __future__ import annotations

import typer


def _predict_tm(seq: str) -> float:
    """PLACEHOLDER Tm (°C). Replace with TemStaPro / DeepSTABp / ThermoMPNN.

    Crude proxy: more charged/polar residues and length correlate with higher Tm in
    thermophiles — enough to exercise the verifier, not to trust."""
    charged = sum(a in "DEKR" for a in seq) / max(len(seq), 1)
    polar = sum(a in "STNQHY" for a in seq) / max(len(seq), 1)
    return min(max(35.0 + 60.0 * charged + 20.0 * polar, 0.0), 110.0)


def main(seq: str = typer.Option(..., help="protein sequence (one-letter)")) -> None:
    print(f"tm={_predict_tm(seq):.1f}")


if __name__ == "__main__":
    typer.run(main)
