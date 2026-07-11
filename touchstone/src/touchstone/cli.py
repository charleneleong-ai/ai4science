"""touchstone CLI — verify a designed metal-binder from the command line.

    touchstone verify design.pdb --metal Ni2+
    touchstone verify design.cif --metal Cu2+ --deep --json

Any agent can call it over Bash; the MCP server (touchstone-mcp) wraps the same engine.
"""

from __future__ import annotations

import json as _json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from .expression import load_signals, score_provider as expression_score_provider
from .geometry.metalhawk import load_predictions, score_provider as metalhawk_score_provider
from .geometry.precedent import metalpdb_precedent_search
from .reward import rank_structures
from .service import verify_structure
from .thermostability import tm_provider

app = typer.Typer(add_completion=False, help="Generator-agnostic verifier for designed metal binders.")
_COLOR = {"trust": "green", "weak": "yellow", "defer": "red"}


@app.callback()
def _root() -> None:
    """Generator-agnostic verifier for designed metal binders."""  # forces `touchstone verify ...`


@app.command()
def verify(
    structure: Path = typer.Argument(..., help="a metal-coordination structure (.pdb / .cif)"),
    metal: str = typer.Option("Ni2+", help="target metal, e.g. Ni2+ / Cu2+ / Co2+"),
    deep: bool = typer.Option(False, "--deep", help="also run the MLIP relaxation + MD (needs a GPU backend)"),
    stress: bool = typer.Option(False, "--stress", help="also test robustness under acidic-leachate / low-pH stress"),
    precedent: bool = typer.Option(False, "--precedent", help="enable the open MetalPDB coordination-motif precedent tier"),
    metalhawk_scores: Path = typer.Option(None, "--metalhawk-scores", help="JSON of precomputed MetalHawk predictions (scripts/metalhawk_score.py) — enables the experimental MetalHawk tier"),
    sequence: str = typer.Option("", "--sequence", help="the design's sequence — enables the expression / thermostability tiers (which key by sequence)"),
    expression_scores: Path = typer.Option(None, "--expression-scores", help="JSON {sequence: {pseudo_perplexity, solubility}} — enables the expression tier (needs --sequence)"),
    thermostability_scores: Path = typer.Option(None, "--thermostability-scores", help="JSON {sequence: Tm_celsius} — enables the thermostability tier (needs --sequence)"),
    json: bool = typer.Option(False, "--json", help="emit JSON (for agents / piping)"),
) -> None:
    """Score a structure's metal site: trust / weak / defer + per-verifier breakdown."""
    metalhawk_scorer = metalhawk_score_provider(load_predictions(metalhawk_scores)) if metalhawk_scores else None
    expression_scorer = expression_score_provider(load_signals(expression_scores)) if expression_scores else None
    tm_predictor = tm_provider(_json.loads(Path(thermostability_scores).read_text())) if thermostability_scores else None
    result = verify_structure(
        structure, metal, deep, stress=stress, sequence=sequence,
        precedent_search=metalpdb_precedent_search if precedent else None,
        metalhawk_scorer=metalhawk_scorer,
        expression_scorer=expression_scorer,
        thermostability_predictor=tm_predictor,
    )
    if json:
        typer.echo(_json.dumps(result, indent=2))
        return

    console = Console()  # created here so CliRunner/redirected stdout is captured
    table = Table(title=f"touchstone · {result['metal']} · {Path(result['structure']).name}")
    for col in ("verifier", "verdict", "score", "reason"):
        table.add_column(col)
    for name, r in result["verifiers"].items():
        if "label" in r:
            table.add_row(name, f"[{_COLOR[r['label']]}]{r['label']}[/]", f"{r['score']:.3f}", r["reason"])
        else:
            table.add_row(name, "[dim]skipped[/]", "—", f"[dim]{r.get('skipped') or r.get('error')}[/]")
    console.print(table)
    c = result["consensus"]
    console.print(
        f"consensus: [{_COLOR[c]}]{c.upper()}[/]  "
        f"(CN {result['coordination_number']}, donors {result['donors'] or '—'})"
    )
    if result.get("stress"):
        conditions = ", ".join(
            f"{cond} [{_COLOR[v['label']]}]{v['label']}[/]" for cond, v in result["stress"].items()
        )
        console.print(f"robustness: {conditions}")
    if result.get("not_run"):
        console.print(f"[dim]not run (needs input): {', '.join(result['not_run'])}[/]")


@app.command()
def rank(
    structures: list[Path] = typer.Argument(..., help="designs to rank (.pdb/.cif; shell-glob expands)"),
    metal: str = typer.Option("Ni2+", help="target metal, e.g. Ni2+ / Cu2+ / Co2+"),
    deep: bool = typer.Option(False, "--deep", help="also run the MLIP relaxation + MD (needs a GPU backend)"),
    precedent: bool = typer.Option(False, "--precedent", help="fold the open MetalPDB precedent tier into the reward"),
    top: int = typer.Option(0, "--top", help="show only the top N (0 = all)"),
    json: bool = typer.Option(False, "--json", help="emit JSON (for agents / piping)"),
) -> None:
    """Rank a batch of designs by verifier reward, best first (the best-of-N selection step)."""
    ranked = rank_structures(
        structures, metal, deep, precedent_search=metalpdb_precedent_search if precedent else None
    )
    if top:
        ranked = ranked[:top]
    if json:
        typer.echo(_json.dumps(ranked, indent=2))
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
            table.add_row(str(i), name, f"{r['reward']:.3f}", f"[{_COLOR[c]}]{c}[/]", str(r["coordination_number"]))
    console.print(table)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
