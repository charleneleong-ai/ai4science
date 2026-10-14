"""Score predicted binder/TNF-alpha complexes on the objectives BoltzGen does not measure.

Everything is derived from per-residue buried surface area (unbound SASA minus bound SASA),
not from atom- or residue-pair counts. Pair counting was tried first and was badly misleading:
it scaled with interpenetration, so clashing poly-alanine junk out-scored real designs, and the
per-residue tallies were not even physically plausible. Buried area is a direct physical
quantity and degrades gracefully.

interface_area    target-side buried area, A^2. A real minibinder interface is ~600-1000.
epitope_focus     share of the target's buried area on the epitope named by --epitope.
cross_reactivity  1 - share of buried area on human/mouse-divergent positions. Objective 2.
ph_switch         interface histidines weighted by their own burial and by standoff from a
                  conserved cationic anchor (R32/K90/K128). Objective 3.

No folding model represents protonation, so ph_switch is a structural proxy, not a free energy
-- and it was refuted: PROPKA on real designs has burial and cation proximity correlating with
the WRONG sign, so the column is kept only to reproduce this campaign's earlier CSVs. Objective
3 is measured by `pka_robustness.py`.

Deployed to the GPU box as `scripts/score_designs.py`, which is the name the nine scripts
there import the chain helpers from; `scripts/score_designs.py` in this repo is a different,
`assess`-based scorer. The box copy is not under version control, so this is its record.

The chain reading, SASA and numbering live in `tnf_binder.structure`; this file holds only the
scoring on top of them. They used to be defined here as well, and the two copies diverged --
this one classified a chain as the target if its identity to human TNF exceeded 0.8, a test
that never fired because mouse TNF scores 0.795.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from tnf_binder.structure import CLASH_DISTANCE as CLASH
from tnf_binder.structure import Complex
from tnf_binder.structure import chain_sequence as chain_seq
from tnf_binder.structure import mature_numbering
from tnf_binder.structure import read_model as load
from tnf_binder.structure import ResidueKey, residue_sasa
from tnf_binder.target import HOTSPOTS, LOBE, POCKET, TnfTarget

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import PH_ANCHORS, divergent_positions  # noqa: E402

# The nine scripts on the box do `from score_designs import ...`; this is that surface. Named
# explicitly so the repo's F401 gate cannot mistake a re-export nothing here calls for a stray
# import and delete the box's dependency.
__all__ = ["CLASH", "chain_seq", "classify_chains", "load", "mature_map", "residue_sasa",
           "low_complexity", "score", "main"]

# Which positions `epitope_focus` is a share of -- a named choice, because the literal that
# used to sit here matched nothing else in the campaign. It was {33,86,87,88,90,91,92,147}, the
# first campaign's generation hotspots, which share no position with the pocket the
# `screen_pocket_*` designs were later generated against (113/115/143/144) and are not the
# conserved TNFR2 site that `submit.py` declares. The three disagree by a lot rather than a
# little: on one submitted design the same structure reads 0.09, 0.18 and 0.56.
EPITOPES: dict[str, tuple[int, ...]] = {
    "conserved": TnfTarget().conserved_epitope,  # what submit.py declares
    "pocket": POCKET,
    "lobe": LOBE,
    "hot": HOTSPOTS,  # the literal this file used to carry, for comparison with its own CSVs
}
DIVERGENT = divergent_positions()
NEAR_ANCHOR = 10.0
HIS_REF_SASA = 180.0
MIN_INTERFACE = 350.0     # below this there is no real interface to speak of
ANCHOR_ATOMS = ("NZ", "NH1", "NH2", "NE")


def classify_chains(model: Model) -> tuple[list[Chain], list[Chain]]:
    """Target chains and binder chains, under the name the other scripts import.

    The rule itself lives in `tnf_binder.structure`. It used to live here too, and the two
    copies diverged: this one also accepted a chain whose identity to human TNF exceeded 0.8,
    which never fired -- mouse TNF scores 0.795 -- so a mouse target classified as the binder.
    The package tests exact containment against both orthologs instead.
    """
    complex_ = Complex.from_model(model)
    return list(complex_.target_chains), list(complex_.binder_chains)


def mature_map(ch: Chain) -> dict[Residue, int]:
    return mature_numbering([ch])


def low_complexity(seq: str) -> bool:
    if not seq:
        return True
    return Counter(seq).most_common(1)[0][1] / len(seq) > 0.5 or len(set(seq)) <= 4


def buried_by_mature_position(targets: Sequence[Chain], t_free: dict[ResidueKey, float],
                              bound: dict[ResidueKey, float],
                              tnum: dict[Residue, int]) -> Counter[int | None]:
    """Target-side buried area per residue, attributed to mature TNF numbering.

    Deliberately not `Complex.buried_by_position`, which drops residues the alignment could not
    place. This keeps them under a `None` key and so counts them in `interface_area`, which
    makes the denominator of every share in the CSV wider than the package's. Swapping in the
    package's version would change every number this file has ever written.
    """
    num_of = {(ch.id, r.id[1]): tnum.get(r) for ch in targets for r in ch if r.id[0] == " "}
    buried: Counter = Counter()
    for key, free in t_free.items():
        delta = free - bound.get(key, free)
        if delta > 0.1:
            buried[num_of.get(key)] += delta
    return buried


def his_term(binders: Sequence[Chain], binder_buried: dict[ResidueKey, float],
             anchors: list[Atom],
             ns: PDB.NeighborSearch) -> tuple[float, list[str], int, int]:
    """(ph_switch, per-histidine detail, interface histidines, clashing atom pairs).

    The clash count rides along because it needs the same walk over the binder's heavy atoms.
    """
    ph, detail, n_his, clashes = 0.0, [], 0, 0
    for ch in binders:
        for r in (x for x in ch if x.id[0] == " "):
            key = (ch.id, r.id[1])
            heavy = [a for a in r if a.element != "H"]
            for a in heavy:
                clashes += sum(1 for _ in ns.search(a.coord, CLASH))
            if r.get_resname() != "HIS" or binder_buried.get(key, 0.0) <= 1.0:
                continue
            n_his += 1
            burial = min(1.0, binder_buried[key] / HIS_REF_SASA)
            d = min((np.linalg.norm(a.coord - b.coord) for a in heavy for b in anchors),
                    default=np.inf)
            prox = max(0.0, 1.0 - d / NEAR_ANCHOR) if np.isfinite(d) else 0.0
            ph += burial * (0.5 + prox)
            detail.append(f"{r.id[1]}:bur={burial:.2f},d={d:.1f}" if np.isfinite(d)
                          else f"{r.id[1]}:bur={burial:.2f},d=NA")
    return ph, detail, n_his, clashes


def flags_for(seq: str, clashes: int, interface_area: float, on_epi: float,
              anchors_found: bool) -> list[str]:
    flags = []
    if clashes:
        flags.append(f"clash({clashes})")
    if low_complexity(seq):
        flags.append("low_complexity")
    if interface_area < MIN_INTERFACE:
        flags.append(f"tiny_interface({interface_area:.0f})")
    if interface_area and on_epi / interface_area < 0.25:
        flags.append("off_epitope")
    if not anchors_found:
        flags.append("anchors_not_found")
    return flags


def score(path: Path, epitope: frozenset[int]) -> dict[str, object]:
    model = load(path)
    targets, binders = classify_chains(model)
    base = {"design": path.stem}
    if not targets or not binders:
        return base | {"flag": "chain_classification_failed"}

    tids = {c.id for c in targets}
    bids = {c.id for c in binders}
    bound = residue_sasa(model, tids | bids)
    t_free = residue_sasa(model, tids)
    b_free = residue_sasa(model, bids)

    tnum: dict[Residue, int] = {}
    for ch in targets:
        tnum.update(mature_map(ch))

    buried_by_pos = buried_by_mature_position(targets, t_free, bound, tnum)
    interface_area = sum(buried_by_pos.values())
    on_epi = sum(v for p, v in buried_by_pos.items() if p in epitope)
    on_div = sum(v for p, v in buried_by_pos.items() if p in DIVERGENT)

    binder_buried = {k: b_free[k] - bound.get(k, b_free[k]) for k in b_free}
    anchors = [a for ch in targets for r in ch if r.id[0] == " "
               and tnum.get(r) in PH_ANCHORS for a in r if a.get_id() in ANCHOR_ATOMS]
    tatoms = [a for ch in targets for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    ph, his_detail, n_his, clashes = his_term(
        binders, binder_buried, anchors, PDB.NeighborSearch(tatoms))

    seq = "".join(chain_seq(ch) for ch in binders)
    flags = flags_for(seq, clashes, interface_area, on_epi, bool(anchors))

    return base | {
        "binder_len": len(seq),
        "interface_area": round(interface_area, 1),
        "epitope_focus": round(on_epi / interface_area, 4) if interface_area else 0.0,
        "divergent_frac": round(on_div / interface_area, 4) if interface_area else 0.0,
        "cross_reactivity": round(1 - on_div / interface_area, 4) if interface_area else 0.0,
        "n_his_interface": n_his,
        "ph_switch": round(ph, 4),
        "clashes": clashes,
        "flag": ",".join(flags),
        "his_detail": ";".join(his_detail),
        "sequence": seq,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--epitope", choices=sorted(EPITOPES), default="conserved",
                    help="which positions epitope_focus is a share of")
    args = ap.parse_args()

    epitope = frozenset(EPITOPES[args.epitope])
    print(f"epitope={args.epitope} {sorted(epitope)}", file=sys.stderr)

    paths = sorted(p for p in args.structures.rglob("*") if p.suffix in (".cif", ".pdb", ".mmcif"))
    if not paths:
        print(f"no structures under {args.structures}")
        return

    # A pipeline output dir holds BOTH the design-step pose and refold_cif/ for the same stem.
    # They are not interchangeable: the design step places grafted scaffolds only approximately
    # (a nanobody measured 524 atom pairs under 2.2 A, interpenetrating at backbone level),
    # while the refold is the model's own physically sensible re-prediction. Scoring a parent dir
    # silently produced two contradictory rows per design, so prefer the refold and say so.
    by_stem: dict[str, Path] = {}
    for q in paths:
        prev = by_stem.get(q.stem)
        if prev is None or ("refold" in q.parts and "refold" not in prev.parts):
            by_stem[q.stem] = q
    if len(by_stem) != len(paths):
        print(f"note: {len(paths)} files collapsed to {len(by_stem)} designs "
              f"(preferred refold_cif where both exist)", file=sys.stderr)
    paths = sorted(by_stem.values())
    rows: list[dict[str, object]] = []
    for p in paths:
        try:
            rows.append(score(p, epitope))
        except Exception as e:  # one bad file must not sink the batch
            rows.append({"design": p.stem, "flag": f"scoring_error:{type(e).__name__}"})
    keys = ["design", "binder_len", "interface_area", "epitope_focus", "cross_reactivity",
            "divergent_frac", "n_his_interface", "ph_switch", "clashes", "flag",
            "his_detail", "sequence"]
    rows.sort(key=lambda r: (bool(r.get("flag")), -float(r.get("ph_switch") or 0),
                             -float(r.get("epitope_focus") or 0)))
    w = csv.DictWriter(sys.stdout, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
    if args.out:
        with args.out.open("w", newline="") as f:
            ww = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            ww.writeheader()
            ww.writerows(rows)
        print(f"\nwrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
