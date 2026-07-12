"""Score a design sequence for the ExpressionVerifier — ESM-2 pseudo-perplexity
(sequence naturalness) + a solubility signal. Runs on a GPU box with `fair-esm`.

Pseudo-perplexity: mask each residue, read the model's log-prob of the true amino
acid, and exponentiate the mean negative log-likelihood. Lower = more natural.

Solubility here is a PLACEHOLDER heuristic (Kyte–Doolittle GRAVY + charge fraction)
— swap in a real predictor (NetSolP / CamSol) before quantitative use; it is the
one signal not yet a learned model.

    conda run -n esm python scripts/expression_score.py --seq MKLV...
"""

from __future__ import annotations

import math

import typer

KD = {  # Kyte–Doolittle hydropathy
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def solubility(seq: str) -> float:
    """PLACEHOLDER: GRAVY (lower → more soluble) blended with charged-residue fraction.
    Returns 0..1. Replace with NetSolP/CamSol."""
    gravy = sum(KD.get(a, 0.0) for a in seq) / max(len(seq), 1)
    charged = sum(a in "DEKR" for a in seq) / max(len(seq), 1)
    return min(max(0.5 - 0.15 * gravy + 0.5 * charged, 0.0), 1.0)


def pseudo_perplexity(seq: str, device: str = "cuda") -> float:
    import esm
    import torch

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    bc = alphabet.get_batch_converter()
    _, _, toks = bc([("x", seq)])
    toks = toks.to(device)
    nll = 0.0
    # one masked forward per residue (serial for clarity); to speed up, stack the L
    # masked copies into one batched forward — same FLOPs, far better GPU utilisation.
    with torch.no_grad():
        for i in range(1, toks.size(1) - 1):  # skip BOS/EOS
            masked = toks.clone()
            masked[0, i] = alphabet.mask_idx
            logits = model(masked)["logits"]
            logp = torch.log_softmax(logits[0, i], dim=-1)
            nll += -logp[toks[0, i]].item()
    return math.exp(nll / max(len(seq), 1))


def main(seq: str = typer.Option(..., help="protein sequence (one-letter)")) -> None:
    ppl = pseudo_perplexity(seq)
    sol = solubility(seq)
    print(f"pseudo_perplexity={ppl:.3f}")
    print(f"solubility={sol:.3f}")


if __name__ == "__main__":
    typer.run(main)
