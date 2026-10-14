"""Assemble the predicted complexes for upload, named as the submission names them.

Proteinbase shows a structure preview, a target, pLDDT and Boltz2 ipSAE per protein, and all
four are blank for this collection. The predicted complexes exist; they are just scattered over
four prediction directories under batch-stamped names that do not match the submitted ones.

Uploads the *complex* rather than the binder alone: ipSAE is an interface score and the target
column wants a target, so a binder on its own would leave both blank again.

Each file is the model_0 sample, which is what every number in the submission's evidence table
was computed on -- the other two diffusion samples are what the robustness analysis used, and
mixing them here would mean the uploaded structure did not match the reported metrics.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

PREDICTIONS = [
    Path("/home/ubuntu/tnf/boltz_sub2_out/boltz_results_boltz_sub2_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_camp_out/boltz_results_boltz_camp_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_fill_out/boltz_results_boltz_fill_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_a38_out/boltz_results_boltz_a38_in/predictions"),
]
COLUMNS = ["name", "file", "source_batch", "source_design", "binder_chain", "binder_length",
           "target_chains", "model0_complex_plddt", "model0_binder_ptm",
           "model0_binder_to_target_iptm", "mean_binder_to_target_iptm",
           "model0_global_iptm_inflated", "n_samples"]


def chain_sequence(chain) -> str:
    return "".join(THREE2ONE.get(r.get_resname(), "X") for r in chain if r.id[0] == " ")


def pair_mean(pair: dict[str, dict[str, float]], binder: str, targets: list[str]) -> float:
    """Symmetrised binder-to-target ipTM, the only number that measures the designed interface."""
    vals = [float(pair[binder][t]) for t in targets] + [float(pair[t][binder]) for t in targets]
    return sum(vals) / len(vals) if vals else 0.0


def main() -> None:
    submission, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    parser = MMCIFParser(QUIET=True)

    index: dict[str, Path] = {}
    for root in PREDICTIONS:
        for directory in sorted(root.iterdir()) if root.exists() else []:
            cif = next(iter(sorted(directory.glob("*_model_0.cif"))), None)
            if cif is None:
                continue
            chains = [(c.id, chain_sequence(c)) for c in parser.get_structure("x", str(cif))[0]]
            chains = [(i, s) for i, s in chains if s]
            if chains:
                index.setdefault(min(chains, key=lambda t: len(t[1]))[1], cif)

    rows = []
    for row in csv.DictReader(submission.open()):
        sequence = row["sequence"].strip().upper()
        cif = index.get(sequence)
        if cif is None:
            print(f"  {row['name']}: NO MATCH — skipped")
            continue
        dest = out / f"{row['name']}.cif"
        shutil.copyfile(cif, dest)

        confidence = json.loads((cif.parent / f"confidence_{cif.stem}.json").read_text())
        # The uploaded file is model_0, but every number the submission ranked on is the MEAN
        # across samples: best-of-n and mean differ by up to 0.35 within one design here, so
        # reporting only the uploaded sample's score would contradict the evidence table.
        samples = sorted(cif.parent.glob("confidence_*.json"))
        model = parser.get_structure("x", str(cif))[0]
        chains = {c.id: chain_sequence(c) for c in model}
        binder = min((i for i, s in chains.items() if s), key=lambda i: len(chains[i]))
        targets = [i for i in confidence["chains_ptm"] if i != "0"]
        b2t = pair_mean(confidence["pair_chains_iptm"], "0", targets)
        per_sample = [pair_mean(json.loads(f.read_text())["pair_chains_iptm"], "0", targets)
                      for f in samples]
        mean_b2t = sum(per_sample) / len(per_sample) if per_sample else 0.0
        rows.append({
            "name": row["name"], "file": dest.name,
            "source_batch": cif.parent.parent.parent.parent.name,
            "source_design": cif.parent.name, "binder_chain": binder,
            "binder_length": len(chains[binder]),
            "target_chains": "".join(i for i in chains if i != binder),
            "model0_complex_plddt": round(confidence["complex_plddt"], 3),
            "model0_binder_ptm": round(confidence["chains_ptm"]["0"], 3),
            "model0_binder_to_target_iptm": round(b2t, 3),
            "mean_binder_to_target_iptm": round(mean_b2t, 3),
            "model0_global_iptm_inflated": round(confidence["iptm"], 3),
            "n_samples": len(per_sample),
        })
        flag = "  <-- model_0 overstates" if b2t - mean_b2t > 0.1 else ""
        print(f"  {row['name']:<9} <- {cif.parent.name:<34} "
              f"model_0 {b2t:.3f}  mean {mean_b2t:.3f}{flag}")

    with (out / "manifest.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{len(rows)} complexes + manifest.csv -> {out}")


if __name__ == "__main__":
    main()
