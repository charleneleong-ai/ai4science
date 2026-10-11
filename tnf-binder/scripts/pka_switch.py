"""Score pH-dependent release using PROPKA-predicted histidine pKa values.

This replaces the structural proxy used earlier (burial x distance-to-cation), which could only
rank candidates. The pH switch has a precise physical condition, so predict it directly:

Objective 3 wants binding at pH 7.4 and none at pH 6.0. A histidine delivers that only if its
pKa sits BETWEEN the two, so the residue is neutral at 7.4 (binding) and cationic at 6.0
(released). pKa below ~6.0 stays neutral at both pHs and never fires; pKa above ~7.4 is
protonated at both and simply blocks binding outright.

This is also why "put a His next to K90" was the wrong instinct: an adjacent cation DEPRESSES
His pKa, pushing it below the trigger. PROPKA computes that electrostatic shift explicitly
instead of guessing at it from distance.

Scoring per histidine:
  ideal    6.3 <= pKa <= 7.2  -- switches cleanly inside the assay window
  usable   6.0 <= pKa <  6.3 or 7.2 < pKa <= 7.6
  dead     outside 6.0-7.6
A design also needs the histidine AT the interface, else protonation changes nothing.
"""
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from score_designs import load, classify_chains, residue_sasa

IDEAL_LO, IDEAL_HI = 6.3, 7.2
USABLE_LO, USABLE_HI = 6.0, 7.6
MIN_INTERFACE_DSASA = 10.0

PKA_LINE = re.compile(r"^\s*(HIS|ASP|GLU|LYS|ARG|CYS|TYR)\s+(\d+)\s+(\S)\s+([-\d.]+)")


def propka_exe() -> str:
    """Resolve propka3 next to the running interpreter; it is not on the subprocess PATH."""
    cand = Path(sys.executable).with_name("propka3")
    return str(cand) if cand.exists() else "propka3"


def run_propka(pdb: Path, workdir: Path) -> dict[tuple[str, int], float]:
    """residue (chain, resnum) -> predicted pKa, for ionizable residues."""
    local = workdir / pdb.name
    local.write_bytes(pdb.read_bytes())
    r = subprocess.run([propka_exe(), local.name], cwd=workdir,
                       capture_output=True, text=True, timeout=1200)
    out = local.with_suffix(".pka")
    if not out.exists():
        print(f"  propka produced no .pka for {pdb.name}: {r.stderr.strip()[:200]}",
              file=sys.stderr)
        return {}
    pkas: dict[tuple[str, int], float] = {}
    in_summary = False
    for line in out.read_text().splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if in_summary:
            if line.strip().startswith("-") or not line.strip():
                if pkas:
                    break
                continue
            m = PKA_LINE.match(line)
            if m:
                _, num, chain, pka = m.groups()
                try:
                    pkas[(chain, int(num))] = float(pka)
                except ValueError:
                    pass
    return pkas


def classify(pka: float) -> str:
    if IDEAL_LO <= pka <= IDEAL_HI:
        return "ideal"
    if USABLE_LO <= pka <= USABLE_HI:
        return "usable"
    return "dead"


def score(pdb: Path, workdir: Path) -> dict:
    model = load(pdb)
    targets, binders = classify_chains(model)
    row = {"design": pdb.stem}
    if not binders:
        return row | {"note": "no binder chain identified"}

    bids = {c.id for c in binders}
    tids = {c.id for c in targets}
    bound = residue_sasa(model, tids | bids) if tids else {}
    b_free = residue_sasa(model, bids)

    pkas = run_propka(pdb, workdir)
    if not pkas:
        return row | {"note": "propka failed"}

    details, n_ideal, n_usable = [], 0, 0
    for ch in binders:
        for r in (x for x in ch if x.id[0] == " "):
            if r.get_resname() != "HIS":
                continue
            key = (ch.id, r.id[1])
            pka = pkas.get(key)
            if pka is None:
                continue
            dsasa = b_free.get(key, 0.0) - bound.get(key, b_free.get(key, 0.0))
            at_interface = dsasa >= MIN_INTERFACE_DSASA
            verdict = classify(pka)
            if at_interface:
                if verdict == "ideal":
                    n_ideal += 1
                elif verdict == "usable":
                    n_usable += 1
            details.append(f"{ch.id}{r.id[1]}:pKa={pka:.2f},dSASA={dsasa:.0f},"
                           f"{verdict}{'' if at_interface else ',not_interface'}")

    return row | {
        "n_his": sum(1 for ch in binders for r in ch if r.get_resname() == "HIS"),
        "n_ideal_switch": n_ideal,
        "n_usable_switch": n_usable,
        "pka_switch_score": round(n_ideal + 0.4 * n_usable, 3),
        "his_detail": ";".join(details),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", type=Path, nargs="+")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    paths: list[Path] = []
    for s in args.structures:
        paths.extend(sorted(p for p in s.rglob("*") if p.suffix == ".pdb")
                     if s.is_dir() else [s])
    if not paths:
        print("no .pdb inputs (propka needs PDB, not mmCIF)")
        return

    rows = []
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        for p in paths:
            rows.append(score(p, wd))
            print(rows[-1])

    if args.out and rows:
        keys = sorted({k for r in rows for k in r})
        with args.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["design"] + [k for k in keys if k != "design"])
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
