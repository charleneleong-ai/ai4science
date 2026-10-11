"""tnf-binder CLI — run the verifier stack on designed TNF-alpha binders.

    tnf-binder epitope                      # the conserved epitope, and what R31->Q costs
    tnf-binder verify design.cif            # every verifier's verdict on one design
    tnf-binder validate submission/         # submission CSVs against Proteinbase's format
    tnf-binder confidence predictions/      # folded designs ranked by binder-to-target ipTM
    tnf-binder wandb-log assets/            # the whole campaign to Weights & Biases

Mirrors touchstone's CLI shape, and shares its Verdict vocabulary, so the two read the same way.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import typer
from rich import box
from rich.console import Console
from rich.table import Table
from touchstone.core import Verdict

from tnf_binder.assess import assess
from tnf_binder.confidence import DesignConfidence, rank_designs
from tnf_binder.submission import csv_problems, row_problem
from tnf_binder.target import LOBE, POCKET, TNFR2_EPITOPE, TnfTarget
from tnf_binder.tracking import CREDIBLE_INTERFACE, FUNNELS, CampaignLogger

VERDICT_COLOR = {"trust": "green", "weak": "yellow", "defer": "red"}
JSON_HELP = "emit JSON (for agents / piping)"

# Re-exported so `tnf_binder.cli` stays the one place a reader looks for the pre-flight rules,
# while `tnf_binder.submission` is the one place they are stated.
__all__ = ["app", "main", "csv_problems", "row_problem"]

app = typer.Typer(add_completion=False)
console = Console()


@app.callback()
def root() -> None:
    """Verifier stack for conditional TNF-alpha binder designs."""
    # Exists only so that the bare command errors instead of silently doing nothing.


def verdict_table(title: str, rows: list[tuple[str, Verdict]]) -> Table:
    table = Table(title=title)
    for column in ("verifier", "verdict", "score", "reason"):
        table.add_column(column)
    for name, verdict in rows:
        table.add_row(name, f"[{VERDICT_COLOR[verdict.label]}]{verdict.label}[/]",
                      f"{verdict.score:.3f}", verdict.reason)
    return table


def interface_bar(value: float, width: int = 16) -> str:
    """A proportional bar for an ipTM, coloured by whether it clears the credibility bar."""
    filled = min(width, round(value * width))
    colour = "green" if value >= CREDIBLE_INTERFACE else "yellow"
    return f"[{colour}]{'█' * filled}[/][grey37]{'╌' * (width - filled)}[/]"


def confidence_table(designs: list[DesignConfidence]) -> Table:
    """Designs ranked on mean binder-to-target ipTM, ruled where they stop being credible."""
    table = Table(
        title="binder-to-target ipTM",
        caption="mean of the symmetrised binder/target pair ipTM across diffusion samples · "
                f"~{CREDIBLE_INTERFACE:.2f} is a credible de novo interface (rule marks the "
                "crossing)\ninflation = global − binder→target, i.e. the target folding "
                "against itself · never rank on global",
        box=box.SIMPLE_HEAD, caption_justify="left",
    )
    table.add_column("#", justify="right", style="grey37")
    table.add_column("design", overflow="fold")  # never ellipsise the identifier
    table.add_column("binder→target", justify="right")
    table.add_column("", justify="left")  # the bar reads as part of the column before it
    for column in ("best-of-n", "global", "inflation", "binder pTM"):
        table.add_column(column, justify="right")

    ruled = False
    for rank, design in enumerate(designs, start=1):
        if not ruled and design.binder_to_target < CREDIBLE_INTERFACE:
            table.add_section()
            ruled = True
        table.add_row(
            str(rank), design.name, f"{design.binder_to_target:.3f}",
            interface_bar(design.binder_to_target),
            f"[grey37]{design.binder_to_target_best:.3f}[/]",
            f"[grey37]{design.global_iptm:.3f}[/]", f"[red]+{design.inflation:.3f}[/]",
            f"{design.binder_ptm:.3f}",
        )
    return table


def epitope_table(target: TnfTarget) -> Table:
    table = Table(title="TNFR2 epitope · 3ALQ contacts within 4.5 A")
    for column in ("position", "human", "mouse", "conserved"):
        table.add_column(column)
    for position in TNFR2_EPITOPE:
        mouse = target.human_to_mouse.get(position)
        mouse_aa = target.mouse.residue(mouse) if mouse else "-"
        conserved = target.human.residue(position) == mouse_aa
        table.add_row(str(position), target.human.residue(position), mouse_aa,
                      "[green]yes[/]" if conserved else "[red]NO[/]")
    return table


@app.command()
def epitope(as_json: bool = typer.Option(False, "--json", help=JSON_HELP)) -> None:
    """Show the TNFR2 epitope, its human/mouse divergence, and the mouse-equivalent hotspots."""
    target = TnfTarget()
    if as_json:
        typer.echo(json.dumps({
            "identity": round(target.identity, 4),
            "epitope_positions": len(TNFR2_EPITOPE),
            "conserved": list(target.conserved_epitope),
            "divergent": sorted(p for p in TNFR2_EPITOPE if p in target.divergent),
            "pocket_human": list(POCKET),
            "pocket_mouse": list(target.mouse_equivalent(POCKET)),
            "lobe_human": list(LOBE),
            "lobe_mouse": list(target.mouse_equivalent(LOBE)),
        }, indent=2))
        return

    console.print(f"human/mouse soluble-domain identity  [bold]{target.identity:.1%}[/]"
                  f"  ({round(target.identity * 157)}/157)")
    console.print(epitope_table(target))
    console.print("\n[bold]R31->Q is the trap[/]: largest SASA in the epitope (193 A^2), beside "
                  "R32 in the 29-33 loop adalimumab and infliximab engage — which is why "
                  "neither is mouse cross-reactive.")


@app.command()
def verify(
    structure: Path = typer.Argument(..., help="a predicted binder/TNF complex (.cif / .pdb)"),
    as_json: bool = typer.Option(False, "--json", help=JSON_HELP),
) -> None:
    """Run every structural verifier on one design: trust / weak / defer, with reasons."""
    assessment = assess(structure)
    if not assessment.is_resolved:
        typer.echo(f"could not separate binder from target in {structure.name}", err=True)
        raise typer.Exit(code=2)
    design, results = assessment.design, list(assessment.verdicts)

    if as_json:
        typer.echo(json.dumps({
            "design": design.name, "length": design.length,
            "verifiers": {n: {"label": v.label, "score": v.score, "reason": v.reason,
                              "metrics": v.metrics} for n, v in results}}, indent=2))
        return

    console.print(verdict_table(f"tnf-binder · {structure.name} · {design.length} aa", results))
    console.print(f"[dim]binder[/] {design.sequence}")


@app.command()
def validate(
    directory: Path = typer.Argument(..., help="directory of submission CSVs"),
) -> None:
    """Check submission CSVs against Proteinbase's stated format before uploading."""
    table = Table(title="submission pre-flight")
    for column in ("file", "rows", "lengths", "verdict"):
        table.add_column(column)

    failed = False
    for path in sorted(directory.glob("*.csv")):
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        problems = csv_problems(rows)
        failed = failed or bool(problems)
        lengths = [len(r["sequence"]) for r in rows] or [0]
        table.add_row(path.name, str(len(rows)), f"{min(lengths)}-{max(lengths)}",
                      f"[red]{'; '.join(problems)}[/]" if problems else "[green]OK[/]")

    console.print(table)
    if failed:
        raise typer.Exit(code=1)


@app.command()
def confidence(
    predictions: Path = typer.Argument(..., help="a Boltz-2 predictions/ directory"),
    as_json: bool = typer.Option(False, "--json", help=JSON_HELP),
) -> None:
    """Rank folded designs by binder-to-target ipTM, showing the global-ipTM inflation."""
    designs = rank_designs(predictions)
    if not designs:
        typer.echo(f"no confidence files under {predictions}", err=True)
        raise typer.Exit(code=2)

    if as_json:
        typer.echo(json.dumps([{
            "design": d.name, "samples": len(d.samples),
            "binder_to_target": round(d.binder_to_target, 3),
            "binder_to_target_best": round(d.binder_to_target_best, 3),
            "global_iptm": round(d.global_iptm, 3),
            "target_internal_iptm": round(d.target_internal, 3),
            "inflation": round(d.inflation, 3), "binder_ptm": round(d.binder_ptm, 3),
        } for d in designs], indent=2))
        return

    console.print(confidence_table(designs))


@app.command("wandb-log")
def wandb_log(
    assets: Path = typer.Argument(..., help="directory holding boltz_camp/ and the pKa scans"),
    project: str = typer.Option("tnf-conditional-binder", "--project"),
    predictions: str = typer.Option("boltz_camp", "--predictions",
                                    help="subdirectory of assets holding the folded designs"),
    funnel: str = typer.Option("campaign", "--funnel",
                               help=f"which denominators this run describes: "
                                    f"{', '.join(FUNNELS)}"),
    entity: str | None = typer.Option(None, "--entity"),
    run_name: str | None = typer.Option(None, "--run-name"),
    notes: str | None = typer.Option(None, "--notes"),
) -> None:
    """Log the funnel, confidence tables, pKa scans and folded complexes to Weights & Biases."""
    if funnel not in FUNNELS:
        typer.echo(f"unknown funnel {funnel!r}; choose from {', '.join(FUNNELS)}", err=True)
        raise typer.Exit(code=2)
    logger = CampaignLogger(assets, predictions_dir=predictions, funnel=FUNNELS[funnel],
                            project=project, entity=entity)
    console.print(f"logging [bold]{len(logger.designs)}[/] folded designs and "
                  f"[bold]{len(logger.structures)}[/] structures to [cyan]{project}[/]")
    try:
        url = logger.log(run_name=run_name, notes=notes)
    except ImportError:
        typer.echo("wandb is not installed: uv sync --extra wandb", err=True)
        raise typer.Exit(code=3) from None
    console.print(f"[green]run[/] {url}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
