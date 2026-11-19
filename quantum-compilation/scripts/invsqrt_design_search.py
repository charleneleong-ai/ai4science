"""Search seed-plus-series invsqrt designs for the lowest modelled Q x C, every input checked.

The model (`quantum_compilation.invsqrt_model`) is calibrated against one measured circuit, and
two of its constants are not pinned down, so the search reports the best variant under each
setting of them rather than one number. Method and findings: docs/invsqrt-diagnosis.md,
"Redesign cost model".

    uv run python scripts/invsqrt_design_search.py <out.csv> --rows <rows.jsonl> --verifiers <dir>
"""

from __future__ import annotations

import csv
import dataclasses
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Annotated, Any

import typer
from rich.progress import track

from quantum_compilation.cli import VerifiersOption
from quantum_compilation.invsqrt import ReferenceValues
from quantum_compilation.invsqrt_model import (
    REFERENCE,
    Calibration,
    CostModel,
    Evaluation,
    SeriesDesign,
    floored,
)
from quantum_compilation.verify import VerificationRun

# What the table holds: (offset, square). Storing the square without the offset never helps.
TABLE_MODES = ((True, True), (True, False), (False, False))
ADDRESS_BITS = range(5, 10)
SEED_BITS = range(7, 12)
TERMS = (2, 3)
FRAC_BITS = range(22, 33)
KICK_FLOORS = range(35, 40)
EXTRA_FRAC_BITS = 2  # also keep this many bits above the smallest feasible precision
OWN_LOOKUP = (CostModel().lookup_toffolis_per_entry, 0.0)  # per entry, fixed: 6c6a734d's synthesis
INIT_LOOKUP = (1.0, -2.0)  # kickmix's init_lookup
DEFAULT_CACHE = Path("~/quantum-compilation-data/runs/invsqrt-model/reference_values.npz")


def cost_models(reference: Evaluation) -> list[CostModel]:
    """The measured constants (adds uncomputed or not), bracketed by 1 and 2 per product, each with
    6c6a734d's synthesised lookup and kickmix's `init_lookup` (2^a - 2 for 2^a entries)."""
    per_product = (
        1.0,
        CostModel.fitted(reference, uncomputed_adds=True).toffolis_per_partial_product,
        CostModel.fitted(reference).toffolis_per_partial_product,
        2.0,
    )
    return [CostModel(p, *lookup) for p, lookup in itertools.product(per_product, (OWN_LOOKUP, INIT_LOOKUP))]


def score(design: SeriesDesign, values: ReferenceValues) -> list[tuple[Evaluation, int]]:
    """The design at every kick floor (one evaluation, re-floored), with each margin."""
    evaluation, phase = design.evaluate(*values.mantissa_exponent)
    return [(evaluation.at_floor(k), values.margin(floored(phase, k))) for k in KICK_FLOORS]


def feasible_variants(config: tuple[Any, ...], cache: Path, verifiers: Path) -> list[tuple[Evaluation, int]]:
    """Feasible floors at the smallest feasible precision and a few bits above it."""
    values = ReferenceValues.load(cache, verifiers)  # cached per worker process
    found: list[tuple[Evaluation, int]] = []
    stop = None
    for f in FRAC_BITS:
        found += [(ev, m) for ev, m in score(SeriesDesign(*config[:3], f, *config[3:]), values) if m > 0]
        if found and stop is None:
            stop = f + EXTRA_FRAC_BITS
        if f == stop:
            break
    return found


def search(cache: Path, verifiers: Path, workers: int) -> list[tuple[Evaluation, int]]:
    configs = [
        (a, s, t, KICK_FLOORS[0], offset, square, tuned)
        for a, s, t, (offset, square), tuned in itertools.product(
            ADDRESS_BITS, SEED_BITS, TERMS, TABLE_MODES, (False, True)
        )
    ]
    found: list[tuple[Evaluation, int]] = []
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(feasible_variants, c, cache, verifiers) for c in configs]
        for future in track(futures, description="configurations"):
            found += future.result()
    return found


def label(design: SeriesDesign) -> str:
    d = design
    return (
        f"A {d.address_bits} S {d.seed_bits} terms {d.terms} F {d.frac_bits} K {d.kick_floor} "
        f"offset {d.store_offset:d} square {d.store_square:d} tuned {d.tuned:d}"
    )


def main(
    out: Annotated[Path, typer.Argument(help="Where to write every feasible variant (.csv).")],
    rows: Annotated[
        Path, typer.Option(help="rows.jsonl of the circuit REFERENCE stands in for (6c6a734d).")
    ],
    verifiers: VerifiersOption,
    cache: Annotated[Path, typer.Option(help="Cached verifier targets.")] = DEFAULT_CACHE,
    workers: Annotated[int, typer.Option()] = 6,
) -> None:
    cache = cache.expanduser()
    measured = VerificationRun.read_rows(rows)["invsqrt"]["detail"]
    qubits, toffolis = measured["num_qubits"], measured["toffoli_count"]
    values = ReferenceValues.load(cache, verifiers)
    reference = next(pair for pair in score(REFERENCE, values) if pair[0].design == REFERENCE)
    found = search(cache, verifiers, workers)
    typer.echo(f"{len(found)} feasible variants; reference {label(REFERENCE)}, measured Q {qubits}")

    default = Calibration.against(reference[0], CostModel.fitted(reference[0]), toffolis, qubits)
    table = [{**dataclasses.asdict(ev.design), **default.cost(ev, m)} for ev, m in found]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    typer.echo("best variant under each cost model (Toffolis per partial product; lookup per entry, fixed):")
    for model in cost_models(reference[0]):
        calibration = Calibration.against(reference[0], model, toffolis, qubits)
        own = CostModel(model.toffolis_per_partial_product, *OWN_LOOKUP)  # 6c6a734d keeps its lookup
        ref_cost = dataclasses.replace(calibration, model=own).cost(*reference)["cost"]
        ev, margin = min(found, key=lambda pair: calibration.cost(*pair)["cost"])
        best = calibration.cost(ev, margin)
        typer.echo(
            f"  {model.toffolis_per_partial_product:.2f}; {model.lookup_toffolis_per_entry:.2f}, "
            f"{model.lookup_fixed_toffolis:+.0f}: "
            f"{label(ev.design)} -> Q {best['qubits']} C {best['C']:.0f}, "
            f"{best['cost'] / ref_cost - 1:+.1%} vs reference"
        )


if __name__ == "__main__":
    typer.run(main)
