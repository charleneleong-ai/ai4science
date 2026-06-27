"""MCP server exposing touchstone's verifier to agents (Claude Code, etc.).

    touchstone-mcp        # runs a stdio MCP server

Exposes one tool, `verify_metal_binder`, backed by the same engine as the CLI.
Requires the optional `mcp` dependency: `pip install touchstone[mcp]`.
"""

from __future__ import annotations

from .service import verify_structure


def main() -> None:
    from mcp.server.fastmcp import FastMCP

    mcp = FastMCP("touchstone")

    @mcp.tool()
    def verify_metal_binder(structure_path: str, metal: str = "Ni2+", deep: bool = False) -> dict:
        """Verify a designed metal-binding site against the touchstone verifier stack.

        Args:
            structure_path: path to the design (.pdb / .cif) containing the metal + protein.
            metal: target metal label, e.g. "Ni2+", "Cu2+", "Co2+".
            deep: also run the MLIP relaxation (needs a GPU backend); default is the
                instant geometry + bond-valence check.

        Returns a dict with per-verifier verdicts and a trust/weak/defer consensus.
        """
        return verify_structure(structure_path, metal, deep)

    mcp.run()


if __name__ == "__main__":
    main()
