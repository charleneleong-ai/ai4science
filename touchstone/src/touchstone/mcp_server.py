"""MCP server exposing touchstone's verifier to agents (Claude Code, etc.).

    touchstone-mcp                         # stdio (local — what Claude Code spawns)
    touchstone-mcp --http --host 0.0.0.0   # streamable-HTTP (host on a GPU box for the deep tiers)

Exposes one tool, `verify_metal_binder`, backed by the same engine as the CLI.
Requires the optional `mcp` dependency: `pip install touchstone[mcp]`.
"""

from __future__ import annotations

import typer

from .service import verify_structure


def serve(
    http: bool = typer.Option(False, "--http", help="serve over streamable-HTTP for remote hosting (default: stdio)"),
    host: str = typer.Option("127.0.0.1", help="bind address for --http (use 0.0.0.0 to expose)"),
    port: int = typer.Option(8000, help="port for --http"),
) -> None:
    from mcp.server.fastmcp import FastMCP  # optional dep — only needed to run the server

    mcp = FastMCP("touchstone", host=host, port=port)

    @mcp.tool()
    def verify_metal_binder(
        structure_path: str, metal: str = "Ni2+", deep: bool = False, stress: bool = False
    ) -> dict:
        """Verify a designed metal-binding site against the touchstone verifier stack.

        Args:
            structure_path: path to the design (.pdb / .cif) containing the metal + protein.
            metal: target metal label, e.g. "Ni2+", "Cu2+", "Co2+".
            deep: also run the MLIP relaxation (needs a GPU backend); default is the
                instant geometry + bond-valence check.
            stress: also re-verify under extreme operating conditions (acidic-leachate bond
                stretch, low-pH donor protonation) → a `stress` map {neutral/leachate/low_pH}.
                Use when the binder must survive a real recovery process, not just stand still.

        Returns a dict with per-verifier verdicts and a trust/weak/defer consensus (plus a
        `stress` robustness map when `stress=True`).
        """
        return verify_structure(structure_path, metal, deep, stress=stress)

    mcp.run(transport="streamable-http" if http else "stdio")


def main() -> None:
    typer.run(serve)


if __name__ == "__main__":
    main()
