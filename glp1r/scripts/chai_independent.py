"""Independent co-folding read, so nothing is ranked by the model that made it.

BoltzGen generates and ranks with Boltz-2 weights, so its own interface score is
a model grading its own family's output. On the TNF-alpha campaign that
self-confirmation was measured rather than assumed: Chai-1 against Boltz-2
correlated at Spearman +0.046 on the quantity that drove the submission ranking
— no relationship at all — while agreeing closely enough on a separate artefact
(the target's internal interfaces inflating a global score) to confirm it was
real. One model's ranking is a reading; two models agreeing is evidence.

chai_lab 0.0.1's `run_inference` computes scores and ranking data but returns
only the output paths, so the inner call is wrapped to capture them. Reported
per chain pair: a global number mixes in the target's own interfaces, which is
the artefact already corrected for on TNF-alpha.

Run inside the Chai container, which needs a writable bind over its downloads
directory:

    apptainer exec --nv -B <work>:/work \\
      -B <writable>:/usr/envs/chai/lib/python3.11/site-packages/downloads \\
      chai.sif python chai_independent.py /work/in /work/out /work/scores.jsonl
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import chai_lab.chai1 as chai1
import numpy as np
import torch

TRUNK_RECYCLES = 3
DIFFUSION_TIMESTEPS = 200
SEED = 0

captured: dict[str, object] = {}
_inner = chai1.run_folding_on_context


def capture(*args: object, **kwargs: object) -> object:
    result = _inner(*args, **kwargs)
    captured["result"] = result
    return result


chai1.run_folding_on_context = capture


def jsonable(value: object) -> object:
    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def already_scored(path: Path) -> set[str]:
    """Jobs with real samples. A job recorded with none is retried, because the
    first TNF run completed all 40 and captured nothing."""
    if not path.exists():
        return set()
    done = set()
    for line in path.open():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("samples"):
            done.add(row["job"])
    return done


def extract(result: object) -> list[dict[str, object]]:
    ranking = result[2] if isinstance(result, tuple) and len(result) >= 3 else None
    if ranking is None:
        return []
    sequence = ranking if isinstance(ranking, (list, tuple)) else [ranking]
    rows: list[dict[str, object]] = []
    for index, rank in enumerate(sequence):
        scores = getattr(rank, "ptm_scores", None)
        row: dict[str, object] = {
            "sample": index,
            "aggregate": jsonable(getattr(rank, "aggregate_score", None)),
        }
        for field in ("complex_ptm", "interface_ptm", "per_chain_pair_iptm", "per_chain_ptm"):
            if scores is not None and hasattr(scores, field):
                row[field] = jsonable(getattr(scores, field))
        rows.append(row)
    return rows


def main() -> None:
    inputs, outputs, scores = (Path(p) for p in sys.argv[1:4])
    outputs.mkdir(parents=True, exist_ok=True)
    done = already_scored(scores)

    for fasta in sorted(inputs.glob("*.fasta")):
        job = fasta.stem
        if job in done:
            continue
        captured.clear()
        started = time.time()
        try:
            chai1.run_inference(
                fasta_file=fasta,
                output_dir=outputs / job,
                num_trunk_recycles=TRUNK_RECYCLES,
                num_diffn_timesteps=DIFFUSION_TIMESTEPS,
                seed=SEED,
                device=torch.device("cuda:0"),
                use_esm_embeddings=True,
            )
        except Exception as exc:  # noqa: BLE001 - one bad design must not stop the scan
            print(f"{job} FAILED {type(exc).__name__}: {str(exc)[:160]}", flush=True)
            continue
        record = {
            "job": job,
            "elapsed_s": round(time.time() - started, 1),
            "samples": extract(captured.get("result")),
        }
        with scores.open("a") as handle:
            handle.write(json.dumps(record) + "\n")
        print(f"{job} done {record['elapsed_s']:.0f}s samples={len(record['samples'])}", flush=True)

    print("CHAI_DONE", flush=True)


if __name__ == "__main__":
    main()
