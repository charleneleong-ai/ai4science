"""Boltz specs for the combinations excluded by a proxy this campaign later disproved.

Six A38H triples hold all three histidines in the pKa window on their design structure, and all
six were dropped without being folded. The reason was not a measurement of them: `helix0_E38H`,
the *single*, put 59% of its interface on human/mouse-divergent residues, and the triples were
assumed to inherit that. Two doubles, `R21H_E67H` and `R21H_R72H`, went the same way on
conserved-epitope fractions of 0.48 and 0.52.

That proxy has since been shown wrong where it mattered most: the `free_1` family carried the
campaign's *highest* conserved-epitope coverage (0.79) and retained only 41-47% of its human
interface when actually folded against the mouse trimer. Cross-reactivity is now measured by
folding against P06804 rather than inferred from which residues the interface touches, and the
eight slot-filling triples already measured that way span 56-98% retained on one chassis -- so
the interface does vary between combinations, and inheritance from a single was never safe.

Eight candidates for seven free slots. Human and mouse complexes written together so both come
from one numbering check.
"""

import sys
from pathlib import Path

HUMAN_TNF = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHV"
             "LLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESG"
             "QVYFGIIAL")
MOUSE_TNF = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVL"
             "LTHTVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQ"
             "VYFGVIAL")
HELIX0 = ("ELEERIEKTAEEVRELVLEGRVAEAVALAVEEICRVGEEVGVPVERRNAVIREFERKLGELAEELEENPEKREE"
          "IAEKVAELVLKTLREIL")

PARENT = {21: "R", 38: "E", 67: "E", 70: "E", 72: "R"}

# the six A38H triples, then the two doubles -- every combination dropped on the proxy
CANDIDATES = [(21, 38, 67), (21, 38, 70), (21, 38, 72),
              (38, 67, 70), (38, 67, 72), (38, 70, 72),
              (21, 67), (21, 72)]


def mutate(positions: tuple[int, ...]) -> str:
    residues = list(HELIX0)
    for position in positions:
        assert residues[position - 1] == PARENT[position], (
            f"position {position}: parent has {residues[position - 1]}, "
            f"expected {PARENT[position]}")
        residues[position - 1] = "H"
    return "".join(residues)


def name(positions: tuple[int, ...]) -> str:
    return "helix0_" + "_".join(f"{PARENT[p]}{p}H" for p in positions)


def spec(binder: str, target: str) -> str:
    chains = [("X", binder)] + [(cid, target) for cid in "ABC"]
    body = "".join(f"  - protein:\n      id: {cid}\n      sequence: {seq}\n      msa: empty\n"
                   for cid, seq in chains)
    return f"version: 1\nsequences:\n{body}"


def main() -> None:
    human_dir, mouse_dir = Path(sys.argv[1]), Path(sys.argv[2])
    for directory in (human_dir, mouse_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for positions in CANDIDATES:
        binder = mutate(positions)
        assert len(binder) == len(HELIX0), f"{name(positions)}: length changed"
        assert binder.count("H") == len(positions), f"{name(positions)}: wrong His count"
        (human_dir / f"A__{name(positions)}.yaml").write_text(spec(binder, HUMAN_TNF))
        (mouse_dir / f"AM__{name(positions)}.yaml").write_text(spec(binder, MOUSE_TNF))
        print(f"{name(positions):28s} {len(binder)} aa, {len(positions)} His")
    print(f"\n{len(CANDIDATES)} human -> {human_dir}\n{len(CANDIDATES)} mouse -> {mouse_dir}")


if __name__ == "__main__":
    main()
