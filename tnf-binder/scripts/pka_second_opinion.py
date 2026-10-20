"""Second pKa opinion from an independent method.

Every pH-switch claim in this campaign rests on PROPKA, which is an empirical
scheme: it assigns shifts from desolvation, hydrogen bonding and Coulomb terms
with fitted coefficients. If its systematic error happens to point the same way
for every histidine at this interface, the in-window rates inherit that error
and the sample spread will not reveal it, because re-running the same method on
more structures measures structural variability, not method bias.

pKAI is a graph neural network trained on experimental and PHMD pKa values. It
shares no machinery with PROPKA, so agreement between them is evidence the
window placement is real and disagreement localises which positions are not
safe to claim.

Reported per position and per sample so the comparison carries the same
denominator as everything else.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import statistics
import sys
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from Bio.PDB import MMCIFParser, PDBIO  # noqa: E402
from pkai.pKAI import pKAI  # noqa: E402

WINDOW = (6.0, 7.6)

# The binder is written as chain X in every spec; the target trimer is A/B/C.
BINDER_CHAIN = "X"


def to_pdb(source: Path, work: Path) -> Path:
    plain = source
    if source.suffix == ".gz":
        plain = work / source.stem
        with gzip.open(source, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    structure = MMCIFParser(QUIET=True).get_structure("x", str(plain))
    pdb = work / f"{plain.stem}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb))
    return pdb


def histidines(pdb: Path) -> dict[int, float]:
    return {
        int(number): float(value)
        for chain, number, name, value in pKAI(str(pdb))
        if name == "HIS" and chain == BINDER_CHAIN
    }


def measure(design_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        for cif in sorted(design_dir.glob("*_model_*.cif*")):
            try:
                pkas = histidines(to_pdb(cif, work))
            except Exception as exc:  # noqa: BLE001 - one bad sample must not stop the scan
                print(f"    {cif.name}: {type(exc).__name__}", flush=True)
                continue
            for position, value in sorted(pkas.items()):
                rows.append({
                    "design": design_dir.name,
                    "sample": cif.stem.rsplit("_", 1)[-1].replace(".cif", ""),
                    "position": position,
                    "pka_pkai": round(value, 2),
                    "in_window": WINDOW[0] <= value <= WINDOW[1],
                })
    return rows


def summarise(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["design"]), int(row["position"])), []).append(row)

    out: list[dict[str, object]] = []
    for (design, position), samples in sorted(grouped.items()):
        values = [float(s["pka_pkai"]) for s in samples]
        inside = sum(1 for s in samples if s["in_window"])
        out.append({
            "design": design,
            "position": position,
            "n_samples": len(samples),
            "pka_mean": round(statistics.mean(values), 2),
            "pka_min": round(min(values), 2),
            "pka_max": round(max(values), 2),
            "spread": round(max(values) - min(values), 2),
            "in_window_samples": f"{inside}/{len(samples)}",
            "in_window_rate": round(inside / len(samples), 2),
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

    print(f"\n{'design':32s} {'His':>4s} {'mean':>6s} {'spread':>7s} {'in-window':>10s}")
    for row in summary:
        print(f"{row['design']:32s} {row['position']:4d} {row['pka_mean']:6.2f}"
              f" {row['spread']:7.2f} {row['in_window_samples']:>10s}", flush=True)
    print(f"\nwrote {out_rows} and {out_summary}")


if __name__ == "__main__":
    main()
