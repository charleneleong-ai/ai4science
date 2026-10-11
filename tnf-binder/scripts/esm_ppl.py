"""ESM-2 masked pseudo-perplexity for the submitted designs.

Lower perplexity means the sequence looks more like the protein language model's idea of a
protein, which is the cheapest available proxy for "will express". It is a proxy and nothing
more: it is uncalibrated against any wet-lab outcome here, and on this set eleven of thirteen
sequences are two to four substitutions apart, so it can separate the two chassis and little
else. Reported as a ranking, never as an expressibility call.

Batched rather than serial -- the L masked copies of one sequence go through as one forward,
which is the same FLOPs at far better GPU utilisation (touchstone's own scorer notes this).
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import esm
import torch

KD = {  # Kyte-Doolittle hydropathy, for the GRAVY term
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHUNK = 24  # masked copies per forward; sized for the ~14 GB the box leaves free


def pseudo_perplexity(model, alphabet, batch_converter, seq: str, device: str) -> float:
    _, _, toks = batch_converter([("x", seq)])
    toks = toks.to(device)
    positions = list(range(1, toks.size(1) - 1))  # skip BOS/EOS
    nll = 0.0
    with torch.no_grad():
        for start in range(0, len(positions), CHUNK):
            block = positions[start:start + CHUNK]
            masked = toks.repeat(len(block), 1)
            for row, pos in enumerate(block):
                masked[row, pos] = alphabet.mask_idx
            logits = model(masked)["logits"]
            for row, pos in enumerate(block):
                logp = torch.log_softmax(logits[row, pos].float(), dim=-1)
                nll += -logp[toks[0, pos]].item()
    return math.exp(nll / len(seq))


def main() -> None:
    rows = list(csv.DictReader(Path(sys.argv[1]).open()))
    out = Path(sys.argv[2])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    if device == "cuda":
        model = model.half()
    batch_converter = alphabet.get_batch_converter()
    print(f"ESM-2 650M on {device}, {len(rows)} sequences", flush=True)

    results = {}
    for row in rows:
        name, seq = row["name"], row["sequence"].strip().upper()
        ppl = pseudo_perplexity(model, alphabet, batch_converter, seq, device)
        gravy = sum(KD.get(a, 0.0) for a in seq) / len(seq)
        results[name] = {"length": len(seq), "pseudo_perplexity": round(ppl, 3),
                         "gravy": round(gravy, 3), "cysteines": seq.count("C")}
        print(f"  {name:9s} len={len(seq):3d} ppl={ppl:8.3f} gravy={gravy:+.2f}", flush=True)

    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
