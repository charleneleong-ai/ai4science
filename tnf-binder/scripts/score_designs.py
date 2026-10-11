"""Score predicted binder/TNF-alpha complexes on the objectives BoltzGen does not measure.

Everything is derived from per-residue buried surface area (unbound SASA minus bound SASA),
not from atom- or residue-pair counts. Pair counting was tried first and was badly misleading:
it scaled with interpenetration, so clashing poly-alanine junk out-scored real designs, and the
per-residue tallies were not even physically plausible. Buried area is a direct physical
quantity and degrades gracefully.

interface_area    target-side buried area, A^2. A real minibinder interface is ~600-1000.
epitope_focus     share of the target's buried area contributed by the intended conserved lobe.
cross_reactivity  1 - share of buried area on human/mouse-divergent positions. Objective 2.
ph_switch         interface histidines weighted by their own burial and by standoff from a
                  conserved cationic anchor (R32/K90/K128). Objective 3.

No folding model represents protonation, so ph_switch is a structural proxy, not a free energy.
"""
import argparse
import copy
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import Align, PDB
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

sys.path.insert(0, str(Path(__file__).parent))
from tnf_common import HUMAN, PH_ANCHORS, divergent_positions

EPITOPE = {86, 87, 88, 90, 91, 92, 33, 147}
DIVERGENT = divergent_positions()
CLASH = 2.2
NEAR_ANCHOR = 10.0
HIS_REF_SASA = 180.0
MIN_INTERFACE = 350.0     # below this there is no real interface to speak of
ANCHOR_ATOMS = ("NZ", "NH1", "NH2", "NE")

ALIGNER = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)


def load(path: Path):
    p = PDB.MMCIFParser(QUIET=True) if path.suffix in (".cif", ".mmcif") else PDB.PDBParser(QUIET=True)
    return p.get_structure(path.stem, path)[0]


def chain_seq(ch) -> str:
    return "".join(THREE2ONE.get(r.get_resname(), "X") for r in ch if r.id[0] == " ")


def identity_to_human(s: str) -> float:
    if not s:
        return 0.0
    a, b = (str(x) for x in ALIGNER.align(HUMAN, s)[0])
    return sum(1 for x, y in zip(a, b) if x == y and x != "-") / max(1, len(s))


def contains_tnf(seq: str) -> bool:
    """True if the chain contains a long exact stretch of human TNF-alpha.

    Substring containment, not fractional identity: BindCraft merges the two target protomers
    into a SINGLE 304-residue chain, so identity/len(chain) is only ~0.5 and an identity
    threshold misclassifies the target as the binder -- which silently attributed target
    histidines to the binder and zeroed every interface dSASA.
    """
    return any(seq[i:i + 25] in HUMAN for i in range(0, max(1, len(seq) - 25), 5))


def classify_chains(model) -> tuple[list, list]:
    target, binder = [], []
    for ch in model:
        s = chain_seq(ch)
        if len(s) < 5:
            continue
        is_target = contains_tnf(s) or identity_to_human(s) > 0.8
        (target if is_target else binder).append(ch)
    return target, binder


def mature_map(ch) -> dict:
    """Target residue object -> mature human resnum, by alignment.

    Needed because BoltzGen renumbers target chains to 1..N, and 1TNF carries Leu where
    UniProt P01375 has Asp143, so neither stored numbering nor substring matching is safe.
    """
    res = [r for r in ch if r.id[0] == " "]
    s = "".join(THREE2ONE.get(r.get_resname(), "X") for r in res)
    a, b = (str(x) for x in ALIGNER.align(HUMAN, s)[0])
    out, hi, si = {}, 0, 0
    for ca, cb in zip(a, b):
        if ca != "-":
            hi += 1
        if cb != "-":
            if ca != "-":
                out[res[si]] = hi
            si += 1
    return out


def residue_sasa(model, keep_ids: set[str]) -> dict[tuple[str, int], float]:
    """Per-residue SASA with only `keep_ids` chains present."""
    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.id not in keep_ids:
            sub.detach_child(ch.id)
    ShrakeRupley().compute(sub, level="R")
    return {(ch.id, r.id[1]): r.sasa for ch in sub for r in ch if r.id[0] == " "}


def low_complexity(seq: str) -> bool:
    if not seq:
        return True
    return Counter(seq).most_common(1)[0][1] / len(seq) > 0.5 or len(set(seq)) <= 4


def quality_flags(seq: str, clashes: int, interface_area: float,
                  on_epitope: float, has_anchors: bool) -> list[str]:
    """Artefacts that disqualify a design outright, rather than merely rank it lower."""
    flags = []
    if clashes:
        flags.append(f"clash({clashes})")
    if low_complexity(seq):
        flags.append("low_complexity")
    if interface_area < MIN_INTERFACE:
        flags.append(f"tiny_interface({interface_area:.0f})")
    if interface_area and on_epitope / interface_area < 0.25:
        flags.append("off_epitope")
    if not has_anchors:
        flags.append("anchors_not_found")
    return flags


def score(path: Path) -> dict:
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

    tnum: dict = {}
    for ch in targets:
        tnum.update(mature_map(ch))
    num_of = {(ch.id, r.id[1]): tnum.get(r) for ch in targets for r in ch if r.id[0] == " "}

    # Buried area per target residue, attributed to mature numbering.
    buried_by_pos: Counter = Counter()
    for key, free in t_free.items():
        d = free - bound.get(key, free)
        if d > 0.1:
            buried_by_pos[num_of.get(key)] += d
    interface_area = sum(buried_by_pos.values())
    on_epi = sum(v for p, v in buried_by_pos.items() if p in EPITOPE)
    on_div = sum(v for p, v in buried_by_pos.items() if p in DIVERGENT)

    # Binder-side burial, used for the histidine term.
    binder_buried = {k: b_free[k] - bound.get(k, b_free[k]) for k in b_free}

    anchors = [a for ch in targets for r in ch if r.id[0] == " "
               and tnum.get(r) in PH_ANCHORS for a in r if a.get_id() in ANCHOR_ATOMS]

    tatoms = [a for ch in targets for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    ns = PDB.NeighborSearch(tatoms)

    ph, his_detail, n_his = 0.0, [], 0
    clashes = 0
    for ch in binders:
        for r in (x for x in ch if x.id[0] == " "):
            key = (ch.id, r.id[1])
            heavy = [a for a in r if a.element != "H"]
            for a in heavy:
                clashes += sum(1 for ta in ns.search(a.coord, CLASH))
            if r.get_resname() != "HIS" or binder_buried.get(key, 0.0) <= 1.0:
                continue
            n_his += 1
            burial = min(1.0, binder_buried[key] / HIS_REF_SASA)
            d = min((np.linalg.norm(a.coord - b.coord) for a in heavy for b in anchors),
                    default=np.inf)
            prox = max(0.0, 1.0 - d / NEAR_ANCHOR) if np.isfinite(d) else 0.0
            ph += burial * (0.5 + prox)
            his_detail.append(f"{r.id[1]}:bur={burial:.2f},d={d:.1f}" if np.isfinite(d)
                              else f"{r.id[1]}:bur={burial:.2f},d=NA")

    seq = "".join(chain_seq(ch) for ch in binders)
    flags = quality_flags(seq, clashes, interface_area, on_epi, bool(anchors))

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
    args = ap.parse_args()

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
    rows = []
    for p in paths:
        try:
            rows.append(score(p))
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
