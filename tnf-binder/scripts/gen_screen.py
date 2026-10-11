"""Emit the config screen for lifting design_to_target_iptm.

Every architecture so far plateaued at d2t_iptm 0.11-0.30 against a ~0.5 bar. Rather than scale
a config that does not work, screen the untested levers on a small budget and scale the winner.

Levers, in expected order of impact:
  epitope      The original hotspots spanned 8 residues, two protomers and 18 A -- far more than
               the 3-5 tight hotspots that are standard practice, and centred on the LEAST
               designable conserved cluster (concavity 81, 1 apolar). The 113/115/143-147 pocket
               is more concave (104) and more hydrophobic (3 apolar), still 100% conserved, and
               still a genuine TNFR2 contact site.
  helix        All-helical binders are markedly easier to design and to re-predict confidently
               than mixed or beta topologies.
  length       A shallow groove may need more binder to grip it.

Indices in the emitted YAML are 1-based CHAIN POSITIONS, read off the structure by
`StoredNumbering` and asserted against the expected residue identities below -- an unguarded
version of this arithmetic already produced a silent five-residue shift once.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.specs import StoredNumbering

app = typer.Typer(add_completion=False)
console = Console()

# Candidate epitopes, mature numbering, with the identities the target file must show there.
EPITOPES: dict[str, tuple[tuple[int, ...], str]] = {
    # concave hydrophobic conserved pocket -- the designability winner
    "pocket": ((113, 115, 143, 144), "PYLF"),
    # the original lobe, narrowed from 8 spread residues to 3 tight ones
    "lobe3": ((88, 90, 92), "QKN"),
    # pocket plus two more conserved contacts, still one protomer
    "pocket5": ((67, 113, 115, 143, 144), "QPYLF"),
}

# name, epitope, length range, helix-biased
CONFIGS: tuple[tuple[str, str, str, bool], ...] = (
    ("pocket_free", "pocket", "60..95", False),
    ("pocket_helix", "pocket", "60..95", True),
    ("pocket5_helix", "pocket5", "70..110", True),
    ("lobe3_free", "lobe3", "60..95", False),
)


def helix_bundle(n: int) -> str:
    """A three-helix bundle pattern of exactly n residues.

    Uniform 'h' is not a designable topology -- a real bundle needs loops to turn between
    helices. Three helices with two short loops is the classic de novo minibinder fold and is
    much easier to re-predict confidently than a mixed or beta topology.
    """
    loop = 4
    helix_total = n - 2 * loop
    outer = helix_total // 3
    middle = helix_total - 2 * outer
    return "h" * outer + "l" * loop + "h" * middle + "l" * loop + "h" * outer


def secondary_structure(length_range: str) -> str:
    """The bundle pattern, sized to the MINIMUM of the length range.

    BoltzGen rejects a secondary_structure string longer than the sampled sequence (its error
    says "shorter than the sequence", but the check is len(spec) > n). Trailing positions
    default to 'u'.
    """
    n_min = int(length_range.split("..")[0])
    pattern = helix_bundle(n_min)
    if len(pattern) != n_min:
        raise ValueError(f"ss pattern {len(pattern)} != {n_min}")
    return pattern


def spec_text(name: str, target: str, epitope: str, binding: str, length_range: str,
              pattern: str | None) -> str:
    ss_line = "" if pattern is None else f'\n      secondary_structure: "{pattern}"'
    ss_desc = "unconstrained" if pattern is None else f"3-helix bundle over {len(pattern)} residues"
    return f"""# TNF-alpha binder screen config: {name}
# epitope={epitope}  length={length_range}  secondary_structure={ss_desc}
#
# Screening to lift design_to_target_iptm off the 0.11-0.30 plateau. Indices are 1-based CHAIN
# POSITIONS; target files keep 1TNF numbering so position = mature - 5, asserted by gen_screen.py.
entities:
  - protein:
      id: X
      sequence: {length_range}{ss_line}
  - file:
      path: {target}
      include:
        - chain: {{id: A}}
        - chain: {{id: B}}
        - chain: {{id: C}}
      binding_types:
        - chain: {{id: A, binding: "{binding}"}}
      structure_groups: "all"
"""


@app.command()
def main(
    structure: Path = typer.Argument(..., help="the trimer target PDB (tnf_trimer.pdb)"),
    out_dir: Path = typer.Argument(..., help="directory to write the screen specs into"),
    target: str = typer.Option("/home/ubuntu/tnf/targets/tnf_trimer.pdb", "--target",
                               help="the target path BoltzGen will read on the design box"),
) -> None:
    """Write one spec per screen config, each with its numbering asserted."""
    numbering = StoredNumbering.from_file(structure)
    out_dir.mkdir(parents=True, exist_ok=True)

    for name, epitope, length_range, helix in CONFIGS:
        positions, expect = EPITOPES[epitope]
        binding = numbering.ranges("A", positions, expect)
        pattern = secondary_structure(length_range) if helix else None
        console.print(f"[bold]{name}[/] epitope={epitope} mature={list(positions)} -> "
                      f"positions {binding} ({expect}) length={length_range} "
                      f"ss={'bundle' if helix else 'free'}")
        (out_dir / f"screen_{name}.yaml").write_text(
            spec_text(name, target, epitope, binding, length_range, pattern))

    console.print(f"wrote [bold]{len(CONFIGS)}[/] screen specs to [cyan]{out_dir}[/]")


if __name__ == "__main__":
    app()
