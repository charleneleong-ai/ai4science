"""`vqcc` -- measure circuits, gate a submission directory, and price the blend."""

from __future__ import annotations

import sys
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import kickmix as km
import typer
from rich.console import Console
from rich.progress import Progress
from rich.table import Table

from quantum_compilation.cost import (
    PROVISIONAL_CALIBRATION,
    TASKS,
    MetricCalibration,
    UnpricedRotationError,
    dfthcsos_cost,
    dfthcsos_rejections,
    dqi_cost,
    invsqrt_cost,
    score,
    task_ratio,
    toffoli_equivalents,
)
from quantum_compilation.resources import Resources
from quantum_compilation.submission import CIRCUIT_SUFFIXES, SubmissionAudit, audit

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()

CalibrationOption = Annotated[
    Path | None,
    typer.Option(
        "--calibration", help="Calibration JSON; omit for the provisional one."
    ),
]


@dataclass(frozen=True)
class VerifierInputs:
    """The quantities only the official verifier can supply.

    Qubit and Toffoli counts come off the circuits, but the dqi failure rate and the
    chemistry task's physics and whole-block-encoding totals do not: the submitted ANGLE
    circuits are a part of what dfthcsos is scored on, not all of it.
    """

    dqi_failures: int | None = None
    dqi_shots: int = 102_400
    e_gap_ha: float | None = None
    lcu_lambda: float | None = None
    hamiltonian_error_ha: float | None = None
    dfthcsos_toffolis: float | None = None
    dfthcsos_qubits: int | None = None

    @property
    def dqi_failure_rate(self) -> float | None:
        return None if self.dqi_failures is None else self.dqi_failures / self.dqi_shots


def load_calibration(path: Path | None) -> MetricCalibration:
    return (
        PROVISIONAL_CALIBRATION if path is None else MetricCalibration.from_json(path)
    )


def warn_if_provisional(calibration: MetricCalibration) -> None:
    if not calibration.measured:
        console.print(f"[yellow]Uncalibrated metric[/] -- {calibration.provenance}")


def circuit_paths(paths: list[Path]) -> list[Path]:
    """Expand directories to the circuit files inside them, keeping the order given."""
    found: list[Path] = []
    for path in paths:
        if path.is_dir():
            found.extend(
                sorted(p for p in path.rglob("*") if p.suffix in CIRCUIT_SUFFIXES)
            )
        else:
            found.append(path)
    if not found:
        raise typer.BadParameter(
            f"no {' or '.join(CIRCUIT_SUFFIXES)} files under {paths}"
        )
    return found


def measured(paths: list[Path], label: str) -> Iterator[tuple[Path, Resources]]:
    """Parse each circuit in turn, reporting progress -- a 1.46 GB circuit takes a while.

    A redrawing bar is unreadable in a nohup log, so only a TTY gets one; anything else
    gets a timestamped line per circuit.
    """
    if sys.stdout.isatty():
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task(label, total=len(paths))
            for path in paths:
                progress.update(task, description=f"{label} {path.name}")
                yield path, Resources.from_file(path)
                progress.advance(task)
        return
    for index, path in enumerate(paths, start=1):
        start = time.monotonic()
        yield path, Resources.from_file(path)
        console.print(
            f"[{time.strftime('%H:%M:%S')}] {label}: {index}/{len(paths)} "
            f"{path.name} ({time.monotonic() - start:.1f}s)"
        )


def equivalents_cell(res: Resources, calibration: MetricCalibration) -> str:
    try:
        return f"{toffoli_equivalents(res, calibration):,.1f}"
    except UnpricedRotationError:
        return "[yellow]unpriced rotations[/]"


def task_costs(
    report: SubmissionAudit,
    by_path: dict[Path, Resources],
    calibration: MetricCalibration,
    inputs: VerifierInputs,
) -> tuple[dict[str, float], dict[str, str]]:
    """Each task's cost, and for the rest why it cannot be computed here."""
    costs: dict[str, float] = {}
    blocked: dict[str, str] = {}

    if report.invsqrt is not None:
        res = by_path[report.invsqrt]
        try:
            costs["invsqrt"] = invsqrt_cost(
                res.num_qubits, toffoli_equivalents(res, calibration)
            )
        except UnpricedRotationError as exc:
            blocked["invsqrt"] = str(exc)
    else:
        blocked["invsqrt"] = "no invsqrt circuit"

    rate = inputs.dqi_failure_rate
    if report.dqi is None:
        blocked["dqi"] = "no dqi circuit"
    elif rate is None:
        blocked["dqi"] = "needs --dqi-failures from the official verifier"
    else:
        res = by_path[report.dqi]
        # T here is the Toffoli count, not Toffoli-equivalents: only invsqrt prices T gates in.
        costs["dqi"] = dqi_cost(res.num_qubits, res.num_toffoli_like, rate)

    chemistry_inputs = (
        inputs.dfthcsos_qubits,
        inputs.dfthcsos_toffolis,
        inputs.e_gap_ha,
        inputs.lcu_lambda,
    )
    if all(value is not None for value in chemistry_inputs):
        qubits, toffolis, e_gap, lcu_lambda = chemistry_inputs
        costs["dfthcsos"] = dfthcsos_cost(qubits, toffolis, e_gap, lcu_lambda)
    else:
        blocked["dfthcsos"] = (
            "needs --dfthcsos-qubits, --dfthcsos-toffolis, --e-gap and --lambda; the "
            "submitted ANGLE circuits are only part of what this task is scored on"
        )
    return costs, blocked


def print_costs(
    costs: dict[str, float], blocked: dict[str, str], calibration: MetricCalibration
) -> None:
    table = Table("task", "cost", "reference", "ratio")
    for task in TASKS:
        if task in costs:
            table.add_row(
                task,
                f"{costs[task]:,.6g}",
                f"{calibration.reference_costs[task]:,.6g}",
                f"{task_ratio(task, costs[task], calibration):.4f}",
            )
        else:
            table.add_row(
                task, "[yellow]—[/]", f"{calibration.reference_costs[task]:,.6g}", "—"
            )
    console.print(table)
    for task, reason in blocked.items():
        console.print(f"[yellow]{task}[/] {reason}")
    if not blocked:
        console.print(f"\n[bold]score {score(costs, calibration):.5f}[/]")


@app.command()
def resources(
    paths: Annotated[list[Path], typer.Argument(help="Circuit files or directories.")],
    calibration_path: CalibrationOption = None,
) -> None:
    """Measured qubits, Toffoli-likes and rotations for each circuit."""
    calibration = load_calibration(calibration_path)
    table = Table("circuit", "qubits", "Toffoli-like", "T", "rotations >T", "C")
    for path, res in measured(circuit_paths(paths), "measuring"):
        table.add_row(
            path.name,
            f"{res.num_qubits:,}",
            f"{res.num_toffoli_like:,}",
            f"{res.num_t_like:,}",
            f"{res.num_fine_rotations:,}",
            equivalents_cell(res, calibration),
        )
    console.print(table)
    warn_if_provisional(calibration)


@app.command()
def check(
    root: Annotated[Path, typer.Argument(help="Submission directory to gate.")],
    calibration_path: CalibrationOption = None,
    dqi_failures: Annotated[
        int | None,
        typer.Option(help="Failing shots the official dqi verifier reported."),
    ] = None,
    dqi_shots: Annotated[
        int, typer.Option(help="Shots the dqi verifier ran.")
    ] = 102_400,
    e_gap: Annotated[
        float | None, typer.Option(help="E_gap in Ha, from the verifier.")
    ] = None,
    lcu_lambda: Annotated[
        float | None,
        typer.Option("--lambda", help="LCU normalization, from the verifier."),
    ] = None,
    hamiltonian_error: Annotated[
        float | None, typer.Option(help="Encoded-Hamiltonian error in Ha.")
    ] = None,
    dfthcsos_toffolis: Annotated[
        float | None, typer.Option(help="Whole-block-encoding Toffoli count.")
    ] = None,
    dfthcsos_qubits: Annotated[
        int | None, typer.Option(help="Whole-block-encoding qubit count.")
    ] = None,
) -> None:
    """Gate a submission directory: layout, ANGLE widths, and every cost its inputs allow."""
    calibration = load_calibration(calibration_path)
    report = audit(root)
    for problem in report.problems:
        console.print(f"[red]layout[/] {problem}")
    if report.ok:
        console.print(
            f"[green]layout ok[/] -- invsqrt, dqi and "
            f"{len(report.angle_circuits)} angle circuit(s) with a factorization"
        )

    paths = [
        p for p in (report.invsqrt, report.dqi, *report.angle_circuits) if p is not None
    ]
    if not paths:
        raise typer.Exit(code=1)
    by_path = dict(measured(paths, "measuring"))

    inputs = VerifierInputs(
        dqi_failures=dqi_failures,
        dqi_shots=dqi_shots,
        e_gap_ha=e_gap,
        lcu_lambda=lcu_lambda,
        hamiltonian_error_ha=hamiltonian_error,
        dfthcsos_toffolis=dfthcsos_toffolis,
        dfthcsos_qubits=dfthcsos_qubits,
    )
    for reason in dfthcsos_rejections(
        angle_circuit_qubits=[by_path[p].num_qubits for p in report.angle_circuits],
        hamiltonian_error_ha=inputs.hamiltonian_error_ha,
        e_gap_ha=inputs.e_gap_ha,
        lcu_lambda=inputs.lcu_lambda,
    ):
        console.print(f"[red]reject[/] {reason}")

    print_costs(*task_costs(report, by_path, calibration, inputs), calibration)
    warn_if_provisional(calibration)


@app.command()
def blend(
    invsqrt: Annotated[float, typer.Option(help="invsqrt cost: Q * C.")],
    dqi: Annotated[float, typer.Option(help="dqi cost: R(eps) * Q^2 * T.")],
    dfthcsos: Annotated[float, typer.Option(help="dfthcsos cost: lambda_eff * T * Q.")],
    calibration_path: CalibrationOption = None,
) -> None:
    """Per-task ratios and the geometric-mean score for three already-measured costs."""
    calibration = load_calibration(calibration_path)
    print_costs({"invsqrt": invsqrt, "dqi": dqi, "dfthcsos": dfthcsos}, {}, calibration)
    warn_if_provisional(calibration)


@app.command()
def calibration(
    write: Annotated[
        Path | None,
        typer.Option(help="Write the active calibration here as an editable JSON."),
    ] = None,
    calibration_path: CalibrationOption = None,
) -> None:
    """Show the active metric calibration, or write it out to fill in with measured numbers."""
    active = load_calibration(calibration_path)
    table = Table("task", "reference cost")
    for task in TASKS:
        table.add_row(task, f"{active.reference_costs[task]:,.6g}")
    console.print(table)
    console.print(
        f"rotation T price: {active.rotation_t_cost or '[yellow]unknown[/]'}\n"
        f"measured: {active.measured}\nprovenance: {active.provenance}"
    )
    if write is not None:
        active.to_json(write)
        console.print(
            f"\nwrote {write} -- replace the costs from `sample_submission/`, which scores 1.0"
        )


@app.command()
def diagram(
    path: Annotated[Path, typer.Argument(help="Circuit file.")],
    out: Annotated[
        Path | None, typer.Option(help="Output file; defaults beside the input.")
    ] = None,
    flame: Annotated[
        bool, typer.Option(help="Write the SVG flame chart instead of HTML.")
    ] = False,
) -> None:
    """Dump kickmix's own HTML diagram or SVG flame chart for a circuit."""
    circuit = km.Circuit.from_file(path)
    content = circuit.flame_chart_svg() if flame else circuit.html_diagram()
    # Release the parsed circuit before encoding the output: on a GB circuit, holding both
    # plus write_text's encoded copy at once is the difference between finishing and swapping.
    del circuit
    destination = out or path.with_suffix(".flame.svg" if flame else ".html")
    destination.write_text(content)
    console.print(f"wrote {destination} ({destination.stat().st_size:,} bytes)")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
