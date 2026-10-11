"""Do the liabilities this campaign reported actually predict expression?

The experiment card lists two risks read off sequence alone and says plainly that neither is
calibrated: eleven of thirteen designs carry an unpaired cysteine, and the two chassis sit at
opposite ends of a charge/hydrophobicity axis. Proteinbase's export carries `expressed` as an
experimental boolean over four thousand designs, which is the first opportunity this campaign
has had to check either claim against a wet-lab outcome rather than asserting it.

Reports base rates and per-feature expression rates, so a feature that does nothing is visible
as a flat line rather than hidden behind a correlation coefficient.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
STANDARD = set(KD)


def truthy(value) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        v = value.strip().lower()
        if v in ("true", "yes", "1"):
            return True
        if v in ("false", "no", "0"):
            return False
    return None


def expressed_call(evaluations: list[dict]) -> bool | None:
    """One verdict per protein: expressed if any replicate says so."""
    calls = [truthy(e.get("value")) for e in evaluations
             if e.get("type") == "experimental" and e.get("metric") == "expressed"]
    calls = [c for c in calls if c is not None]
    return any(calls) if calls else None


def rate(label: str, subset: list[bool]) -> str:
    if not subset:
        return f"  {label:<34}      n=0"
    return f"  {label:<34} {sum(subset) / len(subset):>6.1%}   n={len(subset)}"


def load(path: pathlib.Path) -> list[dict]:
    """Designs carrying both a standard-residue sequence and an expression verdict."""
    data = []
    for r in csv.DictReader(path.open(encoding="utf-8-sig")):
        seq = (r["sequence"] or "").strip().upper()
        if not seq or set(seq) - STANDARD:
            continue
        try:
            call = expressed_call(json.loads(r["evaluations"] or "[]"))
        except Exception:
            continue
        if call is not None:
            data.append({"seq": seq, "expressed": call, "cys": seq.count("C"),
                         "gravy": sum(KD[a] for a in seq) / len(seq), "len": len(seq)})
    return data


def bands(data: list[dict], key: str, edges, label) -> None:
    """One expression rate per band, so a feature that does nothing reads as a flat column."""
    for lo, hi in edges:
        sub = [d["expressed"] for d in data if lo <= d[key] < hi]
        if sub:
            print(rate(label(lo, hi), sub))


def main() -> None:
    data = load(pathlib.Path(sys.argv[1]))
    print(f"{len(data)} designs with a sequence and an expression verdict")
    print(rate("BASE RATE", [d["expressed"] for d in data]))

    print("\ncysteine count — the campaign's headline liability:")
    for n in range(0, 5):
        sub = [d["expressed"] for d in data if d["cys"] == n]
        if sub:
            print(rate(f"{n} cysteine{'s' if n != 1 else ''}", sub))
    odd = [d["expressed"] for d in data if d["cys"] % 2 == 1]
    even = [d["expressed"] for d in data if d["cys"] % 2 == 0 and d["cys"] > 0]
    none = [d["expressed"] for d in data if d["cys"] == 0]
    print(rate("  -> ODD (unpairable)", odd))
    print(rate("  -> EVEN and non-zero", even))
    print(rate("  -> none", none))

    print("\nhydrophobicity (GRAVY) — the aggregation claim:")
    bands(data, "gravy", ((-2.0, -0.4), (-0.4, -0.2), (-0.2, 0.0), (0.0, 0.2), (0.2, 2.0)),
          lambda lo, hi: f"{lo:+.1f} to {hi:+.1f}")

    print("\nlength — the feature that actually moves:")
    bands(data, "len", ((0, 60), (60, 80), (80, 100), (100, 150), (150, 10000)),
          lambda lo, hi: f"{lo}-{hi} aa")

    print("\nthe tnf-binder designs sit at:")
    print("  helix_0 (11 of 13):  1 cysteine (odd), GRAVY -0.37, 91 aa")
    print("  free_1  (2 of 13):   0 cysteines,      GRAVY +0.26, 74 aa")


if __name__ == "__main__":
    main()
