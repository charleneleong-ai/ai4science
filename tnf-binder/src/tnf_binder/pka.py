"""Running PROPKA, which is how the pH switch stopped being a matter of opinion.

`verifiers.pka_switch` records why this is measured rather than reasoned about: three
structural mechanisms were argued for and PROPKA refuted each on real designs. This module is
the mechanics of getting the number.

Two things bite, both of them silent:

  propka3 is not on the subprocess PATH when the caller is a uv- or conda-managed interpreter,
  so it is resolved next to `sys.executable` first.

  propka needs a COMPLETE side chain. A residue renamed to HIS that keeps only backbone + CB
  has no imidazole (ND1/CD2/CE1/NE2); propka omits it from the summary without complaining and
  every pKa comes back missing. Build the side chain with `rosetta.HisMutator` before asking.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

PKA_LINE = re.compile(r"^\s*(HIS|ASP|GLU|LYS|ARG|CYS|TYR)\s+(\d+)\s+(\S)\s+([-\d.]+)")
PROPKA_TIMEOUT = 1200


def propka_executable() -> str:
    """Resolve propka3 next to the running interpreter; it is not on the subprocess PATH."""
    candidate = Path(sys.executable).with_name("propka3")
    return str(candidate) if candidate.exists() else "propka3"


def parse_summary(text: str) -> dict[tuple[str, int], float]:
    """(chain, residue number) -> pKa, from the SUMMARY block of a .pka file."""
    pkas: dict[tuple[str, int], float] = {}
    in_summary = False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        if line.strip().startswith("-") or not line.strip():
            if pkas:  # the rule after the table, rather than the one before it
                break
            continue
        if match := PKA_LINE.match(line):
            _residue, number, chain, pka = match.groups()
            pkas[(chain, int(number))] = float(pka)
    return pkas


def run_propka(pdb: Path, workdir: Path) -> dict[tuple[str, int], float]:
    """Predicted pKa per ionizable residue, or nothing if propka could not score the file.

    propka writes its output beside its input, so the structure is copied into `workdir` first
    rather than littering the caller's directory.
    """
    local = workdir / pdb.name
    local.write_bytes(pdb.read_bytes())
    result = subprocess.run([propka_executable(), local.name], cwd=workdir,
                            capture_output=True, text=True, timeout=PROPKA_TIMEOUT, check=False)
    output = local.with_suffix(".pka")
    if not output.exists():
        print(f"  propka produced no .pka for {pdb.name}: {result.stderr.strip()[:200]}",
              file=sys.stderr)
        return {}
    return parse_summary(output.read_text())
