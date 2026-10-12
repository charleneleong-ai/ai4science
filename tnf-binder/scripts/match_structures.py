"""Find the folded complex for each submitted design, by sequence rather than by name.

The structures are spread over four prediction directories, each stamped with its own batch
marker, and the design names in them do not all match the names in the submission -- the
free_1 variant is `A73H` in the campaign's chain numbering and `T73H` in the evidence table.
Matching on the binder sequence instead of on a name means a mismatch cannot be introduced by
a naming convention, which is the only way an upload could silently attach the wrong structure
to a design.

Prints one line per submitted design: the directory holding it, and whether the binder chain
in that structure is character-for-character the sequence that was submitted.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

PREDICTIONS = [
    Path("/home/ubuntu/tnf/boltz_sub2_out/boltz_results_boltz_sub2_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_camp_out/boltz_results_boltz_camp_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_fill_out/boltz_results_boltz_fill_in/predictions"),
    Path("/home/ubuntu/tnf/boltz_a38_out/boltz_results_boltz_a38_in/predictions"),
]


def chain_sequence(chain) -> str:
    return "".join(THREE2ONE.get(r.get_resname(), "X") for r in chain if r.id[0] == " ")


def binder_of(cif: Path) -> tuple[str, str] | None:
    """(chain id, sequence) of the shortest chain -- the binder against a 157-aa target."""
    parser = MMCIFParser(QUIET=True)
    try:
        model = parser.get_structure(cif.stem, str(cif))[0]
    except Exception:
        return None
    chains = [(c.id, chain_sequence(c)) for c in model]
    chains = [(i, s) for i, s in chains if s]
    return min(chains, key=lambda t: len(t[1])) if chains else None


def main() -> None:
    index: dict[str, tuple[Path, str]] = {}
    for root in PREDICTIONS:
        if not root.exists():
            continue
        for directory in sorted(root.iterdir()):
            if not directory.is_dir():
                continue
            cif = next(iter(sorted(directory.glob("*_model_0.cif"))), None)
            if cif is None:
                continue
            found = binder_of(cif)
            if found:
                index.setdefault(found[1], (directory, found[0]))
    print(f"indexed {len(index)} distinct binder sequences across "
          f"{sum(r.exists() for r in PREDICTIONS)} prediction directories\n")

    rows = list(csv.DictReader(Path(sys.argv[1]).open()))
    matched = 0
    for row in rows:
        sequence = row["sequence"].strip().upper()
        hit = index.get(sequence)
        if hit:
            matched += 1
            directory, chain = hit
            print(f"  {row['name']:<9} chain {chain}  {directory.parent.parent.parent.name}/"
                  f"{directory.name}")
        else:
            near = [(sum(a != b for a, b in zip(s, sequence)), d.name)
                    for s, (d, _c) in index.items() if len(s) == len(sequence)]
            near.sort()
            hint = f" (nearest: {near[0][1]}, {near[0][0]} substitutions)" if near else ""
            print(f"  {row['name']:<9} NO EXACT MATCH{hint}")
    print(f"\n{matched} of {len(rows)} submitted designs have a folded complex "
          f"whose binder sequence matches character-for-character")


if __name__ == "__main__":
    main()
