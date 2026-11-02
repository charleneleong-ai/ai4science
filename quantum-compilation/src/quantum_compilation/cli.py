"""`vqcc` -- gate a submission, run the official verifiers, and price the blend.

The division of labour is the thing to keep straight. `check` is instant and claims only
what filenames and a circuit's width can support. `verify` runs the official verifiers and
is the only command that produces a real cost: no task's cost is computable without them,
invsqrt included, because the scored Toffoli count is a simulated per-shot mean and the
price of a rotation depends on the accuracy margin the simulation leaves over.
"""

from __future__ import annotations

import time
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Annotated

import kickmix as km
import typer
import yaml
from rich.console import Console
from rich.progress import Progress
from rich.table import Table

from quantum_compilation import verify as verifier
from quantum_compilation.chemistry import ChemistryMetrics, ChemistryScreen
from quantum_compilation.cost import (
    DEFAULT_CALIBRATION,
    ESCAPE_SCORE,
    MAX_BROT,
    MIN_BROT,
    TASKS,
    MetricCalibration,
    ShapeResources,
    blend_factor,
    cost_cut_for_blend,
    dfthcsos_rejections,
    lambda_eff_floor,
    marginal_blend_per_unit,
    require_task,
    score,
    task_ratio,
)
from quantum_compilation.submission import (
    CHEMISTRY_DIR,
    CIRCUIT_SUFFIXES,
    Resources,
    SubmissionAudit,
)

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()

CalibrationOption = Annotated[
    Path | None,
    typer.Option("--calibration", help="Calibration JSON; omit for the built-in one."),
]
VerifiersOption = Annotated[
    Path,
    typer.Option(
        "--verifiers",
        envvar="VQCC_VERIFIERS",
        help="The extracted `verifiers/` directory from the competition download.",
    ),
]


def load_calibration(path: Path | None) -> MetricCalibration:
    return DEFAULT_CALIBRATION if path is None else MetricCalibration.from_json(path)


def warn_if_unmeasured(calibration: MetricCalibration) -> None:
    if not calibration.measured:
        console.print(f"[yellow]Uncalibrated metric[/] -- {calibration.provenance}")


def fmt_cost(cost: float | None) -> str:
    return f"{cost:,.8g}" if cost is not None else "[yellow]—[/]"


def circuit_paths(paths: list[Path]) -> list[Path]:
    """Expand directories to the circuit files inside them, keeping the order given."""
    found: list[Path] = []
    for path in paths:
        if path.is_dir():
            found.extend(sorted(p for p in path.rglob("*") if p.suffix in CIRCUIT_SUFFIXES))
        else:
            found.append(path)
    if not found:
        raise typer.BadParameter(f"no {' or '.join(CIRCUIT_SUFFIXES)} files under {paths}")
    return found


def measured(paths: Sequence[Path], label: str) -> Iterator[tuple[Path, Resources]]:
    """Parse each circuit in turn, reporting progress -- a 1.46 GB circuit takes a while.

    A redrawing bar is unreadable in a nohup log, so only a terminal gets one; anything else
    gets a timestamped line per circuit.
    """
    if console.is_terminal:
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


def print_costs(costs: dict[str, float], calibration: MetricCalibration) -> None:
    table = Table("task", "cost", "reference", "ratio")
    for task in TASKS:
        cost = costs.get(task)
        table.add_row(
            task,
            "[yellow]—[/]" if cost is None else fmt_cost(cost),
            fmt_cost(calibration.reference_costs[task]),
            "—" if cost is None else f"{task_ratio(task, cost, calibration):.4f}",
        )
    console.print(table)
    if len(costs) == len(TASKS):
        console.print(f"\n[bold]score {score(costs, calibration):.5f}[/]")
    else:
        console.print(f"[yellow]no score[/] -- needs all of {TASKS}, have {sorted(costs)}")


def share_table(
    unit: str, whole: str, total: float, rows: Sequence[tuple[str, float, str]]
) -> Table:
    """A term/value/share-of-total table, closed by the bolded total itself.

    Shared by the `T` and `Q` tables, which were otherwise the same loop twice. The two
    float-valued share tables above stay as they are: routing them through here would cost a
    value-format and an optional-column parameter, which is more than the repetition.
    """
    table = Table("term", unit, f"share of {whole}", "bought with")
    for name, value, note in rows:
        table.add_row(name, f"{value:,.0f}", f"{value / total:.2%}", note)
    table.add_row(f"[bold]{whole}_total[/]", f"[bold]{total:,.0f}[/]", "", "")
    return table


def print_chem_metrics(metrics: ChemistryMetrics) -> None:
    """The screened metrics, then the terms behind `Lambda` and the gap."""
    table = Table("metric", "value", "reads")
    for name, value, reads in (
        ("error_1", f"{metrics.error_1:.9f}", "one-body residual"),
        ("error_2", f"{metrics.error_2:.9f}", "two-body Frobenius residual"),
        ("E_gap", f"{metrics.e_gap:.9f}", "E_upper - E_SOS"),
        ("Lambda", f"{metrics.one_norm:.6f}", "LCU one-norm"),
        ("lambda_eff", f"{metrics.lambda_eff:.6f}", "sqrt(E_gap * D)"),
        ("D = 2*Lambda - E_gap", f"{metrics.gap_product:.6f}", "what the score turns on"),
        ("0.5 * sum(W^2)", f"{metrics.half_sum_w2:.6f}", "raises E_gap, lowers D"),
        ("c0", f"{metrics.c0:.6f}", "fixed by the reference Hamiltonian"),
    ):
        table.add_row(name, value, reads)
    console.print(table)

    split = Table("term of Lambda", "value", "share")
    for name, value in (
        ("lambda_h2", metrics.lambda_h2),
        ("positive eigen mass", metrics.positive_mass),
        ("negative eigen mass", metrics.negative_mass),
    ):
        split.add_row(name, f"{value:.6f}", f"{value / metrics.one_norm:.2%}")
    console.print(split)


def print_chem_headroom(metrics: ChemistryMetrics) -> None:
    """How much of `D` a factorization search can actually move.

    `D = 2P + 2 * two_body_excess - c0` with both leading terms non-negative, so `-c0` is a
    floor no refit can cross. Printed on every screen because it is the number that decides
    whether a candidate chemistry direction is worth any compute at all.
    """
    table = Table("term of D", "value", "share")
    for name, value in (
        ("movable: 2P + 2*excess", metrics.movable_mass),
        ("fixed: -c0", -metrics.c0),
    ):
        table.add_row(name, f"{value:.6f}", f"{value / metrics.gap_product:.2%}")
    console.print(table)
    console.print(
        f"lambda_eff floor [bold]{lambda_eff_floor(metrics.c0):.6f}[/] against "
        f"{metrics.lambda_eff:.6f} -- the whole factorization search is worth "
        f"[bold]{metrics.lambda_eff_headroom:.3%}[/]"
    )


def print_shape_resources(res: ShapeResources, blend: float, target: float) -> None:
    """Where `T` and `Q` go, then what a point of blend costs in each.

    Printed together on purpose: the shares say which of the two is worth attacking, and the
    prices say whether the attack can reach a stated target at all.
    """
    shape, free = "the shape (R, B, C)", "nothing -- submitted circuits"
    console.print(
        share_table(
            "Toffolis",
            "T",
            res.toffoli,
            [
                ("ANGLE circuits", res.angle_toffoli, free),
                ("inner state prep", res.inner_toffoli, shape),
                ("Givens rotations", res.givens_toffoli, "brot"),
                ("outer state prep", res.outer_toffoli, shape),
            ],
        )
    )
    console.print(
        share_table(
            "qubits",
            "Q",
            res.qubits,
            [
                ("shared workspace", res.heuristic_temporary_qubits, shape),
                ("persistent state prep", res.heuristic_persistent_qubits, shape),
                ("angle-word registers", res.angle_word_qubits, "brot"),
            ],
        )
    )

    cut = cost_cut_for_blend(blend, target)
    console.print(
        f"\none Toffoli [bold]{marginal_blend_per_unit(blend, res.toffoli):+.6f}[/] blend, "
        f"one qubit [bold]{marginal_blend_per_unit(blend, res.qubits):+.6f}[/]"
    )
    console.print(
        f"{blend:.5f} -> {target:g} needs [bold]{cut:.2%}[/] off the chemistry cost: "
        f"{res.toffoli * cut:,.0f} Toffolis, or {res.qubits * cut:,.1f} qubits, or both"
    )
    reach = res.angle_cut_for(cut)
    verdict = "reachable" if reach <= 1 else "[red]out of reach[/]"
    console.print(
        f"from the free share alone: ANGLE loaders [bold]{reach:.2%}[/] cheaper -- {verdict}"
    )


def agreement_cell(result: verifier.Verification) -> str:
    """Whether `cost` reproduced the verifier's score, as three distinct states."""
    match result.model_agrees:
        case True:
            return f"[green]{fmt_cost(result.model_cost)} ✓[/]"
        case False:
            return f"[red]{fmt_cost(result.model_cost)} ✗[/]"
        case _:
            return "[dim]not evaluated[/]"


def selected_tasks(tasks: list[str] | None, whole_submission: bool) -> tuple[str, ...]:
    """The tasks to run, rejecting a bad name or a filter that cannot do the asked-for job."""
    selected = tuple(tasks) if tasks else TASKS
    for task in selected:
        try:
            require_task(task)
        except KeyError as exc:
            raise typer.BadParameter(str(exc)) from exc
    if whole_submission and len(selected) != len(TASKS):
        # Caught before the run rather than after: a calibration needs all three costs, and
        # finding that out at the end would waste the whole simulation.
        raise typer.BadParameter(
            f"--write-calibration needs every task; --task selected {sorted(selected)}"
        )
    return selected


def print_results_table(results: Sequence[verifier.Verification]) -> None:
    table = Table("task", "status", "cost", "model cost")
    for result in results:
        table.add_row(
            result.task,
            "[green]PASS[/]" if result.passed else f"[red]{result.status}[/]",
            fmt_cost(result.score),
            agreement_cell(result),
        )
    console.print(table)


def write_calibration_file(
    run: verifier.VerificationRun,
    destination: Path,
    root: Path,
    salt: str,
    seed: int | None,
) -> None:
    """Write a measured calibration, warning if it is not the known reference submission."""
    if differing := DEFAULT_CALIBRATION.disagreeing_hashes(run.hashes):
        console.print(
            f"[yellow]not the known reference submission[/] -- {sorted(differing)} hash "
            "differently from the circuits the built-in calibration was measured on. "
            "Calibrating against these makes them score 1.0 by construction."
        )
    provenance = (
        f"measured: official verifiers over {root} at salt {salt!r}"
        f"{'' if seed is None else f' seed {seed}'} on {time.strftime('%Y-%m-%d')}"
    )
    run.calibration(provenance).to_json(destination)
    console.print(f"[green]wrote {destination}[/] -- pass it as --calibration")


@app.command()
def resources(
    paths: Annotated[list[Path], typer.Argument(help="Circuit files or directories.")],
) -> None:
    """Static gate counts per circuit -- a triage bound, not a cost.

    Toffoli-like here is the largest number any shot could execute. The scored figure is the
    simulated mean, which measurement-based uncomputation puts below it, so this bounds a
    candidate's cost from above: useful for rejecting one without simulating it, never for
    ranking two that both survive.
    """
    table = Table("circuit", "qubits", "Toffoli-like (max)", "T", "rotations >T")
    for path, res in measured(circuit_paths(paths), "measuring"):
        table.add_row(
            path.name,
            f"{res.num_qubits:,}",
            f"{res.num_toffoli_like:,}",
            f"{res.num_t_like:,}",
            f"{res.num_fine_rotations:,}",
        )
    console.print(table)


@app.command()
def check(
    root: Annotated[Path, typer.Argument(help="Submission directory to gate.")],
) -> None:
    """Gate a submission's layout, and the one rejection knowable without simulating.

    Everything here is cheap and certain. For a cost, run `vqcc verify`.
    """
    report = SubmissionAudit.from_root(root)
    for problem in report.problems:
        console.print(f"[red]layout[/] {problem}")
    if report.ok:
        console.print(
            f"[green]layout ok[/] -- invsqrt, dqi and "
            f"{len(report.angle_circuits)} angle circuit(s) with a factorization"
        )
    if not report.circuits:
        raise typer.Exit(code=1)

    if report.angle_circuits:
        # Only the width is wanted, so read it off the parsed circuit rather than paying for
        # a full operation_statistics() pass per file.
        angle_qubits = [km.Circuit.from_file(path).num_qubits for path in report.angle_circuits]
        reasons = dfthcsos_rejections(angle_circuit_qubits=angle_qubits)
        for reason in reasons:
            console.print(f"[red]reject[/] {reason}")
        if not reasons:
            console.print(
                f"[green]widths ok[/] -- widest ANGLE circuit {max(angle_qubits)} qubits"
            )
    raise typer.Exit(code=0 if report.ok else 1)


@app.command()
def verify(
    root: Annotated[Path, typer.Argument(help="Submission directory to score.")],
    verifiers: VerifiersOption,
    calibration_path: CalibrationOption = None,
    tasks: Annotated[
        list[str] | None, typer.Option("--task", help="Limit to these tasks.")
    ] = None,
    salt: Annotated[
        str, typer.Option(help="Salt merged into the candidate's hash to seed the simulator.")
    ] = verifier.DEFAULT_SALT,
    seed: Annotated[
        int | None,
        typer.Option(
            help="Explicit simulator seed. Pass one to compare candidates: without it each "
            "candidate draws its own seed from its own bytes."
        ),
    ] = None,
    timeout: Annotated[
        float | None, typer.Option(help="Seconds per task before giving up.")
    ] = None,
    write_calibration: Annotated[
        Path | None,
        typer.Option(help="With `sample_submission/`: write the measured calibration here."),
    ] = None,
    rows: Annotated[
        Path | None,
        typer.Option(help="Write one JSON row per verdict here, updated as each lands."),
    ] = None,
) -> None:
    """Run the official verifiers and price the result -- the real score, free and unlimited.

    This is the dev set: the leaderboard runs this same code and rations measurements to
    five a day, so every candidate comes through here and only a winner costs a submission.
    Pass `--seed` when comparing candidates; see the `verify` module for why a salt alone
    does not hold the draw fixed.
    """
    calibration = load_calibration(calibration_path)
    selected = selected_tasks(tasks, whole_submission=write_calibration is not None)

    report = SubmissionAudit.from_root(root)
    for problem in report.problems:
        console.print(f"[yellow]layout[/] {problem}")
    if seed is None:
        console.print(
            "[dim]no --seed: each task is seeded from its own circuit's hash, which is the "
            "published default but not comparable across candidates[/]"
        )

    # The run rewrites `rows` as each verdict lands, so an interrupt late in a ~12-minute
    # simulation keeps everything already measured.
    run = verifier.VerificationRun(rows_path=rows)
    for result in run.measure(verifiers, report, selected, salt=salt, seed=seed, timeout=timeout):
        console.print(f"[dim]{result.task}: {result.status} in {result.seconds:.1f}s[/]")
        if not result.passed:
            console.print(f"[red]{result.task}[/] {result.error}")

    print_results_table(run.results)
    if rows is not None:
        console.print(f"wrote {rows}")

    print_costs(run.costs, calibration)
    if write_calibration is None:
        warn_if_unmeasured(calibration)
    else:
        write_calibration_file(run, write_calibration, root, salt, seed)


@app.command()
def blend(
    invsqrt: Annotated[float, typer.Option(help="invsqrt cost: Q * C.")],
    dqi: Annotated[float, typer.Option(help="dqi cost: R(eps) * Q^2 * T.")],
    dfthcsos: Annotated[float, typer.Option(help="dfthcsos cost: lambda_eff * T * Q.")],
    calibration_path: CalibrationOption = None,
) -> None:
    """Per-task ratios and the geometric-mean score for three already-measured costs."""
    calibration = load_calibration(calibration_path)
    print_costs({"invsqrt": invsqrt, "dqi": dqi, "dfthcsos": dfthcsos}, calibration)
    warn_if_unmeasured(calibration)


@app.command()
def chem(
    factorization: Annotated[Path, typer.Argument(help="A dfthcsos `submission.npz`.")],
    verifiers: VerifiersOption,
    toffolis: Annotated[
        float | None, typer.Option(help="T from a verifier run, to turn lambda_eff into a cost.")
    ] = None,
    qubits: Annotated[int | None, typer.Option(help="Q from a verifier run.")] = None,
) -> None:
    """Screen a chemistry factorization's constraints without building its circuits.

    Fast enough to search with -- under a second a candidate rather than a minute -- because
    it reads the submitted angle words instead of recovering them from a simulation. That is
    also its limit: the verifier scores the words it recovers, and it cannot know T or Q, so
    this reports a cost only when given them.
    """
    metrics = ChemistryScreen(verifiers).evaluate(factorization)
    print_chem_metrics(metrics)
    print_chem_headroom(metrics)
    for reason in metrics.rejected_for:
        console.print(f"[red]reject[/] {reason}")
    if not metrics.rejected_for:
        console.print(
            "[green]screened constraints ok[/] -- ANGLE widths need the circuits, "
            "see `vqcc check`"
        )
    if toffolis is not None and qubits is not None:
        cost = fmt_cost(metrics.cost(toffolis, qubits))
        console.print(f"cost with T={toffolis:,.0f} Q={qubits}: {cost}")


@app.command()
def shape(
    parts: Annotated[
        Path, typer.Argument(help="A submission's parts record, e.g. docs/experiments/*-parts.yaml.")
    ],
    verifiers: VerifiersOption,
    blend: Annotated[float, typer.Option(help="The blend this submission scores.")],
    target: Annotated[float, typer.Option(help="The blend to price the gap to.")] = ESCAPE_SCORE,
    sweep_brot: Annotated[
        bool, typer.Option("--sweep-brot", help="Also cost every brot the model admits.")
    ] = False,
) -> None:
    """Split `T` and `Q` into the part a factorization buys and the part it compiles for free.

    Reads the shape, `brot` and the ANGLE Toffoli mean from the record of what was submitted,
    and derives the word-register count `J` from its `Q_total`. Only the heuristic half is
    modelled, so `--sweep-brot` costs nothing.
    """
    part = yaml.safe_load(parts.read_text())["parts"][CHEMISTRY_DIR]
    screen = ChemistryScreen(verifiers)
    here = screen.resources_from_record(part)
    print_shape_resources(here, blend, target)
    if not sweep_brot:
        return

    # lambda_eff and the ANGLE mean are both held fixed -- the first drops out of a ratio of
    # costs, the second is a property of the loaders rather than of the model -- so this ranks
    # the heuristic half alone, with no chemistry and no circuits.
    shape_args = (int(part["R"]), int(part["B"]), int(part["C"]))
    reference = here.toffoli * here.qubits
    table = Table("brot", "T", "Q", "workspace", "blend", "delta")
    for candidate in range(MIN_BROT, MAX_BROT + 1):
        try:
            res = screen.resources(
                *shape_args,
                candidate,
                angle_toffoli=here.angle_toffoli,
                angle_word_registers=here.angle_word_qubits // int(part["brot"]),
            )
        except ValueError as exc:
            table.add_row(str(candidate), f"[yellow]{exc}[/]", "", "", "", "")
            continue
        moved = blend * blend_factor(res.toffoli * res.qubits, reference)
        table.add_row(
            str(candidate),
            f"{res.toffoli:,.0f}",
            f"{res.qubits:,}",
            f"{res.heuristic_temporary_qubits:,}",
            f"{moved:.5f}",
            f"{moved - blend:+.4f}",
        )
    console.print(table)


@app.command()
def calibration(
    write: Annotated[
        Path | None, typer.Option(help="Write the active calibration here as JSON.")
    ] = None,
    calibration_path: CalibrationOption = None,
) -> None:
    """Show the active metric calibration -- each task's reference cost and where it came from."""
    active = load_calibration(calibration_path)
    table = Table("task", "reference cost", "measured on")
    for task in TASKS:
        digest = active.reference_sha256.get(task, "")
        table.add_row(
            task,
            fmt_cost(active.reference_costs[task]),
            f"{digest[:12]}..." if digest else "—",
        )
    console.print(table)
    console.print(f"measured: {active.measured}\nprovenance: {active.provenance}")
    if write is not None:
        active.to_json(write)
        console.print(
            f"\nwrote {write} -- or let `vqcc verify sample_submission/ "
            "--write-calibration` measure it"
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
