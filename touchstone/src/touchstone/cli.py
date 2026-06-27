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

from .service import verify_structure

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
    json: bool = typer.Option(False, "--json", help="emit JSON (for agents / piping)"),
) -> None:
    """Score a structure's metal site: trust / weak / defer + per-verifier breakdown."""
    result = verify_structure(structure, metal, deep)
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
    if result.get("not_run"):
        console.print(f"[dim]not run (needs input): {', '.join(result['not_run'])}[/]")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
