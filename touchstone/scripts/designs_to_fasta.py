"""Extract protein-chain sequences from design structures → a FASTA for TemStaPro.

TemStaPro judges global fold stability, so it needs the *full* sequence — but touchstone's
parser keeps only the coordination site. This reads each design's protein chain with gemmi
and writes one FASTA record per file (id = filename stem), ready to hand to
`scripts/thermostability_score.py --fasta`. The metal/ligand het groups are dropped (only
the longest polymer chain is kept). Runs in any env with gemmi (the A100 design envs):

    conda run -n bg python scripts/designs_to_fasta.py --glob 'boltzgen_out/**/*.cif' --out designs.fasta
"""

from __future__ import annotations

from pathlib import Path

import gemmi
import typer


def _sequence(path: Path) -> str:
    """The design's protein sequence: the longest polymer chain, one-letter, gaps removed."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    seqs = [c.get_polymer().make_one_letter_sequence() for c in st[0] if len(c.get_polymer())]
    return max(seqs, key=len).replace("-", "").upper() if seqs else ""


def main(
    glob: str = typer.Option(..., help="glob for design structures, e.g. 'designs/*.pdb'"),
    out: Path = typer.Option(..., help="output FASTA path"),
) -> None:
    records = []
    for p in sorted(Path().glob(glob)):
        if seq := _sequence(p):
            records.append((p.stem, seq))
        else:
            print(f"skip {p.name}: no protein chain")
    out.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    print(f"wrote {out} ({len(records)} sequences)")


if __name__ == "__main__":
    typer.run(main)
