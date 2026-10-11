"""Write His-mutant design files that BoltzGen's own folding step can re-score.

Re-using BoltzGen's folding step (rather than a fresh Boltz co-fold) means the mutants are
judged by the same iptm / interaction-pae metrics as the parents, so the comparison is
apples-to-apples. The folding step re-predicts coordinates from the design, so a mutated
residue only needs backbone plus CB to be structurally sane; the remaining His side-chain
atoms are rebuilt during folding.

Emits <design_dir>/<name>.cif plus the parent's .npz metadata alongside, which is the pair
layout FromGeneratedDataModule expects (suffix .cif, suffix_metadata .npz).
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN

# Atoms retained when mutating a residue to histidine. Folding rebuilds the imidazole.
BACKBONE = {"N", "CA", "C", "O", "CB", "OXT"}


def is_tnf(seq: str) -> bool:
    return any(seq[i:i + 25] in HUMAN for i in range(0, max(1, len(seq) - 25), 5))


def binder_chain(st: gemmi.Structure) -> gemmi.Chain | None:
    for ch in st[0]:
        seq = gemmi.one_letter_code([r.name for r in ch])
        if seq and len(seq) < 260 and not is_tnf(seq):
            return ch
    return None


def mutate_to_his(src: Path, positions: list[int], dest: Path) -> bool:
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    ch = binder_chain(st)
    if ch is None:
        return False
    res = [r for r in ch]
    for pos in positions:
        if not 1 <= pos <= len(res):
            return False
        r = res[pos - 1]
        for atom in [a.name for a in r]:
            if atom not in BACKBONE:
                for i, a in enumerate(r):
                    if a.name == atom:
                        del r[i]
                        break
        r.name = "HIS"
    st.setup_entities()
    dest.parent.mkdir(parents=True, exist_ok=True)
    st.make_mmcif_document().write_file(str(dest))
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("variants_csv", type=Path)
    ap.add_argument("parents_dir", type=Path, help="dir holding the parent .cif/.npz pairs")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--limit", type=int, default=60)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.variants_csv.open()))[:args.limit]
    made = 0
    for r in rows:
        parent, name = r["parent"], r["design"]
        positions = [int(x) for x in r["his_positions"].split()]
        src = next((p for p in args.parents_dir.rglob(f"{parent}.cif")), None)
        if src is None:
            print(f"  no parent cif for {parent}")
            continue
        dest = args.out_dir / f"{name}.cif"
        if not mutate_to_his(src, positions, dest):
            print(f"  failed to mutate {name}")
            continue
        npz = src.with_suffix(".npz")
        if npz.exists():
            shutil.copy2(npz, args.out_dir / f"{name}.npz")
        made += 1

    print(f"wrote {made}/{len(rows)} His mutants -> {args.out_dir}")


if __name__ == "__main__":
    main()
