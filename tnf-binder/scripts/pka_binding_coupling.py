"""Is each histidine's titration actually coupled to binding?

The submission gate asks whether a histidine's predicted pKa falls in the
6.0-7.6 window in the bound complex. That is necessary but not sufficient. A
histidine that titrates in the window and titrates *identically* whether the
binder is bound or free contributes nothing to pH-dependent binding: the
protonation equilibrium is the same on both sides, so it cancels.

What governs the pH dependence of Kd is the linkage between the two states. If
the free binder's pKa is higher than the complex's, protonation is favoured in
the free state, so dropping the pH drives dissociation -- the switch the
challenge asks for. If the complex's pKa is higher, protonation stabilises the
complex and binding gets *stronger* at low pH, which is the wrong direction. If
they match, binding is pH-independent at that site.

So the quantity of interest is

    delta = pKa(complex) - pKa(free binder)

and a working switch needs delta meaningfully negative. This was never measured;
only the bound-state pKa was. Run over every diffusion sample so the result
carries the same denominator as the rest of the campaign.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import statistics
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "/home/ubuntu/tnf/scripts")

from Bio.PDB import PDBIO, Select  # noqa: E402
from pka_switch import run_propka  # noqa: E402
from score_designs import classify_chains, load  # noqa: E402

# A histidine needs its pKa to drop by at least this much on binding before the
# coupling is called real rather than numerical noise. PROPKA reproduces to 0.01
# on repeat runs of the same file, but across diffusion samples of the same
# design the spread is far larger, so the threshold is set on the sample spread
# rather than on solver precision. Declared before the run.
COUPLING_FLOOR = 0.3


class OneChain(Select):
    """Write a single chain, so PROPKA sees the binder without its target."""

    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return chain.id == self.chain_id


def unpack(cif: Path, work: Path) -> Path:
    """Structures were gzipped to reclaim disk; PROPKA needs them on disk plain."""
    if cif.suffix != ".gz":
        return cif
    out = work / cif.stem
    with gzip.open(cif, "rb") as src, out.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def pkas_for(model, chain_id: str, work: Path, stem: str, apo: bool) -> dict[int, float]:
    io = PDBIO()
    io.set_structure(model)
    path = work / f"{stem}_{'apo' if apo else 'complex'}.pdb"
    io.save(str(path), OneChain(chain_id) if apo else Select())
    raw = run_propka(path, work)
    return {position: value for (chain, position), value in raw.items() if chain == chain_id}


def measure(design_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        for cif in sorted(design_dir.glob("*_model_*.cif*")):
            plain = unpack(cif, work)
            model = load(plain)
            _target, binders = classify_chains(model)
            if not binders:
                continue
            chain = binders[0]
            positions = [r.id[1] for r in chain if r.get_resname() == "HIS"]
            if not positions:
                continue
            stem = plain.stem
            complexed = pkas_for(model, chain.id, work, stem, apo=False)
            free = pkas_for(model, chain.id, work, stem, apo=True)
            for position in positions:
                bound, unbound = complexed.get(position), free.get(position)
                if bound is None or unbound is None:
                    continue
                rows.append({
                    "design": design_dir.name,
                    "sample": stem.rsplit("_", 1)[-1],
                    "position": position,
                    "pka_complex": round(bound, 2),
                    "pka_apo": round(unbound, 2),
                    "delta": round(bound - unbound, 2),
                })
    return rows


def summarise(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int], list[float]] = {}
    for row in rows:
        grouped.setdefault((str(row["design"]), int(row["position"])), []).append(float(row["delta"]))

    out: list[dict[str, object]] = []
    for (design, position), deltas in sorted(grouped.items()):
        coupled = sum(1 for d in deltas if d <= -COUPLING_FLOOR)
        out.append({
            "design": design,
            "position": position,
            "n_samples": len(deltas),
            "delta_mean": round(statistics.mean(deltas), 2),
            "delta_median": round(statistics.median(deltas), 2),
            "coupled_samples": f"{coupled}/{len(deltas)}",
            "coupled_rate": round(coupled / len(deltas), 2),
            "direction": "release at low pH" if statistics.mean(deltas) <= -COUPLING_FLOOR
            else "wrong direction" if statistics.mean(deltas) >= COUPLING_FLOOR
            else "uncoupled",
        })
    return out


def main() -> None:
    root = Path(sys.argv[1])
    out_rows, out_summary = Path(sys.argv[2]), Path(sys.argv[3])

    all_rows: list[dict[str, object]] = []
    for design_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        rows = measure(design_dir)
        all_rows.extend(rows)
        print(f"{design_dir.name:32s} {len(rows):4d} his-sample measurements", flush=True)

    with out_rows.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)

    summary = summarise(all_rows)
    with out_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    print(f"\n{'design':32s} {'His':>4s} {'dmean':>7s} {'coupled':>9s}  direction")
    for row in summary:
        print(f"{row['design']:32s} {row['position']:4d} {row['delta_mean']:7.2f}"
              f" {row['coupled_samples']:>9s}  {row['direction']}", flush=True)
    print(f"\nwrote {out_rows} and {out_summary}")


if __name__ == "__main__":
    main()
