"""touchstone CLI — verify a designed metal-binder from the command line.

    touchstone verify design.pdb --metal Ni2+
    touchstone verify design.cif --metal Cu2+ --deep --json

Any agent can call it over Bash; the MCP server (touchstone-mcp) wraps the same engine.
"""

from __future__ import annotations

import json as jsonlib
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from .calibration import BINDS, EXPRESSED, bands, calibrate as calibration, load_export
from .expression import load_signals, score_provider as expression_score_provider
from .geometry.metalhawk import load_predictions, score_provider as metalhawk_score_provider
from .reward import rank_structures
from .service import verify_structure
from .thermostability import tm_provider

# The gates worth checking first on a new target: each is a one-line sequence rule that some
# project has argued for from first principles, and all of them are testable in one call.
GATES = {
    "no cysteine": lambda d: "C" not in d.sequence,
    "at most one cysteine": lambda d: d.sequence.count("C") <= 1,
    "even cysteine count": lambda d: d.sequence.count("C") % 2 == 0,
    "under 100 residues": lambda d: len(d.sequence) < 100,
    "under 150 residues": lambda d: len(d.sequence) < 150,
    "hydrophilic (GRAVY < 0)": lambda d: sum(KD_HYDROPATHY.get(a, 0.0) for a in d.sequence) < 0,
}
LENGTH_EDGES = (0, 60, 80, 100, 150, 10_000)
KD_HYDROPATHY = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
                 "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
                 "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

app = typer.Typer(add_completion=False, help="Generator-agnostic verifier for designed metal binders.")
COLOR = {"trust": "green", "weak": "yellow", "defer": "red"}


@app.callback()
def root() -> None:
    """Generator-agnostic verifier for designed metal binders."""  # forces `touchstone verify ...`


@app.command()
def verify(
    structure: Path = typer.Argument(..., help="a metal-coordination structure (.pdb / .cif)"),
    metal: str = typer.Option("Ni2+", help="target metal, e.g. Ni2+ / Cu2+ / Co2+"),
    cutoff: float = typer.Option(2.8, help="metal–donor distance cutoff in Å. Raise it for a long axial donor — a type-1 Cu Met thioether sits at ~2.82 Å and is excluded at the default, turning N2S2 into N2S1"),
    deep: bool = typer.Option(False, "--deep", help="also run the MLIP relaxation + MD (needs a GPU backend)"),
    stress: bool = typer.Option(False, "--stress", help="also test robustness under acidic-leachate / low-pH stress"),
    precedent: bool = typer.Option(True, "--precedent/--no-precedent", help="open MetalPDB coordination-motif precedent tier (on by default)"),
    metalhawk_scores: Path = typer.Option(None, "--metalhawk-scores", help="JSON of precomputed MetalHawk predictions (scripts/metalhawk_score.py) — enables the experimental MetalHawk tier"),
    sequence: str = typer.Option("", "--sequence", help="the design's sequence — enables the composition / expression / thermostability tiers (which key by sequence)"),
    composition: bool = typer.Option(False, "--composition", help="model-free amino-acid composition gate — rejects the degenerate poly-Ala/Gly sequences that interface-confidence models score highly (needs --sequence)"),
    expression_scores: Path = typer.Option(None, "--expression-scores", help="JSON {sequence: {pseudo_perplexity, solubility}} — enables the expression tier (needs --sequence)"),
    thermostability_scores: Path = typer.Option(None, "--thermostability-scores", help="JSON {sequence: Tm_celsius} — enables the thermostability tier (needs --sequence)"),
    selectivity: str = typer.Option("", "--selectivity", help="comma-separated competing metals for the selectivity tiers, e.g. Ni2+,Cu2+,Co2+ — motif enrichment (MetalPDB occupancy) runs on CPU; the MLIP metal-swap ΔE tier is added by --deep (and currently defers: no backbone passes the Irving-Williams gate)"),
    json: bool = typer.Option(False, "--json", help="emit JSON (for agents / piping)"),
) -> None:
    """Score a structure's metal site: trust / weak / defer + per-verifier breakdown."""
    metalhawk_scorer = metalhawk_score_provider(load_predictions(metalhawk_scores)) if metalhawk_scores else None
    expression_scorer = expression_score_provider(load_signals(expression_scores)) if expression_scores else None
    tm_predictor = tm_provider(jsonlib.loads(Path(thermostability_scores).read_text())) if thermostability_scores else None
    result = verify_structure(
        structure, metal, deep, cutoff=cutoff, stress=stress, sequence=sequence,
        precedent=precedent, composition=composition,
        metalhawk_scorer=metalhawk_scorer,
        expression_scorer=expression_scorer,
        thermostability_predictor=tm_predictor,
        selectivity_metals=selectivity.split(",") if selectivity else None,
    )
    if json:
        typer.echo(jsonlib.dumps(result, indent=2))
        return

    console = Console()  # created here so CliRunner/redirected stdout is captured
    table = Table(title=f"touchstone · {result['metal']} · {Path(result['structure']).name}")
    for col in ("verifier", "verdict", "score", "reason"):
        table.add_column(col)
    for name, r in result["verifiers"].items():
        if "label" in r:
            table.add_row(name, f"[{COLOR[r['label']]}]{r['label']}[/]", f"{r['score']:.3f}", r["reason"])
        else:
            table.add_row(name, "[dim]skipped[/]", "—", f"[dim]{r.get('skipped') or r.get('error')}[/]")
    console.print(table)
    c = result["consensus"]
    console.print(
        f"consensus: [{COLOR[c]}]{c.upper()}[/]  "
        f"(CN {result['coordination_number']}, donors {result['donors'] or '—'})"
    )
    if result.get("stress"):
        conditions = ", ".join(
            f"{cond} [{COLOR[v['label']]}]{v['label']}[/]" for cond, v in result["stress"].items()
        )
        console.print(f"robustness: {conditions}")
    if result.get("not_run"):
        console.print(f"[dim]not run (needs input): {', '.join(result['not_run'])}[/]")


@app.command()
def rank(
    structures: list[Path] = typer.Argument(..., help="designs to rank (.pdb/.cif; shell-glob expands)"),
    metal: str = typer.Option("Ni2+", help="target metal, e.g. Ni2+ / Cu2+ / Co2+"),
    cutoff: float = typer.Option(2.8, help="metal–donor distance cutoff in Å — must match the cutoff the pool is being compared against"),
    deep: bool = typer.Option(False, "--deep", help="also run the MLIP relaxation + MD (needs a GPU backend)"),
    precedent: bool = typer.Option(True, "--precedent/--no-precedent", help="fold the open MetalPDB precedent tier into the reward (on by default)"),
    top: int = typer.Option(0, "--top", help="show only the top N (0 = all)"),
    json: bool = typer.Option(False, "--json", help="emit JSON (for agents / piping)"),
) -> None:
    """Rank a batch of designs by verifier reward, best first (the best-of-N selection step)."""
    ranked = rank_structures(
        structures, metal, deep, precedent=precedent, cutoff=cutoff
    )
    if top:
        ranked = ranked[:top]
    if json:
        typer.echo(jsonlib.dumps(ranked, indent=2))
        return

    console = Console()
    table = Table(title=f"touchstone rank · {metal} · {len(ranked)} designs")
    for col in ("#", "design", "reward", "consensus", "CN"):
        table.add_column(col)
    for i, r in enumerate(ranked, 1):
        name = Path(r["structure"]).name
        if "error" in r:
            table.add_row(str(i), name, "0.000", "[dim]error[/]", "—")
        else:
            c = r["consensus"]
            table.add_row(str(i), name, f"{r['reward']:.3f}", f"[{COLOR[c]}]{c}[/]", str(r["coordination_number"]))
    console.print(table)


def main() -> None:
    app()


@app.command()
def calibrate(
    export: Path = typer.Argument(..., help="a Proteinbase CSV export of designs with outcomes"),
    target: str = typer.Option(None, "--target", help="one target slug; base rates differ between targets"),
    outcome: str = typer.Option(EXPRESSED, "--outcome", help=f"{EXPRESSED} or {BINDS}"),
) -> None:
    """Show what the sequence-level gates would have done to designs with measured outcomes.

    Run this before choosing a threshold for a new target, not after. Every number a verifier
    in this package gates on was a guess until there was a control set to check it against.
    """
    console = Console()  # created here so CliRunner/redirected stdout is captured
    designs = [d for d in load_export(export, target=target) if d.measured(outcome) is not None]
    if not designs:
        typer.echo(f"no design in {export} has a measured {outcome}", err=True)
        raise typer.Exit(code=1)
    console.print(f"[bold]{len(designs)}[/] designs with a measured {outcome}"
                  + (f" against {target}" if target else ""))
    if target is None:
        console.print("[yellow]pooled over every target — base rates differ by two orders of "
                      "magnitude between them, and a gate pooled this way can invert its own "
                      "sign. Pass --target before believing anything below.[/]")

    table = Table(title=f"what each gate would have done ({outcome})")
    for column in ("gate", "verdict", "keeps", "drops", "separation", "z"):
        table.add_column(column)
    for label, gate in GATES.items():
        c = calibration(designs, gate, outcome=outcome)
        colour = {"discriminates": "green", "inverted": "red"}.get(c.verdict, "yellow")
        table.add_row(label, f"[{colour}]{c.verdict}[/]",
                      f"{c.kept} at {c.kept_rate:.1%}", f"{c.dropped} at {c.dropped_rate:.1%}",
                      f"{c.separation:+.1%}", f"{c.z:.1f}")
    console.print(table)

    console.print(f"\nlength, before any threshold is chosen (base rate "
                  f"{calibration(designs, lambda _d: True, outcome=outcome).base_rate:.1%}):")
    for band in bands(designs, lambda d: len(d.sequence), LENGTH_EDGES, outcome=outcome):
        console.print(f"  {band}")


if __name__ == "__main__":
    main()
