"""Read the co-fold screen and ask whether each design needs the peptide.

Ranking tier 2 is peptide-dependence. The test is the same binder against the complex and
against the peptide-free receptor: a design that reads the groove scores alike in both, and one
that reads the peptide should lose confidence when the peptide is gone. The quantity reported is
that drop, and the bar is the drop measured for the 7OW6 receptor, which is known to need the
peptide.

The gate is the drop in whole-complex ipTM. `pair_chains_iptm` is reported alongside it and is
deliberately **not** gated: on the control it gives 0.894 for the receptor against beta-2
microglobulin, two chains that make no contact at all, because a pair value rises with
confidence in the whole assembly rather than with contact between the pair. The drop does not
inherit that defect, because it compares two runs containing the same chains -- but note this
reasoning covers the drop only. There is no per-condition difference available for the pair
value itself, since the apo condition has no peptide to pair against, so the pair column is one
uncalibrated number per design and is here for audit rather than for ranking.

Chain order comes from `cofold_specs`, which wrote the specs, rather than from re-reading them:
the order is ours by construction, and an earlier version recovered it by regexing the YAML
this pipeline had just generated.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import typer
from common import rounded, write_rows, wrote
from cofold_specs import (
    APO,
    COMPLEX,
    CONTROL,
    PEPTIDE_CHAIN,
    WILDTYPE,
    binder_ids,
    chain_ids,
    conditions_from,
    parse_run_name,
)

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class Dependence:
    design: str
    complex_iptm: float
    apo_iptm: float
    wildtype_iptm: float | None
    complex_peptide_pair: float | None
    complex_plddt: float

    @property
    def drop(self) -> float:
        return self.complex_iptm - self.apo_iptm

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "iptm_complex_g12d": rounded(self.complex_iptm),
            "iptm_apo_hla": rounded(self.apo_iptm),
            "peptide_dependence_drop": rounded(self.drop),
            "iptm_complex_wt": rounded(self.wildtype_iptm),
            "design_peptide_pair_iptm": rounded(self.complex_peptide_pair),
            "complex_plddt": rounded(self.complex_plddt),
        }


def peptide_pair(
    confidence: dict[str, object], order: tuple[str, ...], binder: tuple[str, ...]
) -> float | None:
    """Mean pair ipTM between the binder chains and the peptide chain, by position in the spec.

    Both directions of each pair are averaged. Boltz's matrix is asymmetric -- the sibling
    campaign measured the two directions differing by up to 0.08 -- and neither is privileged.
    """
    pairs = confidence.get("pair_chains_iptm")
    if not isinstance(pairs, dict) or PEPTIDE_CHAIN not in order:
        return None
    peptide = order.index(PEPTIDE_CHAIN)

    def cell(row: int, column: int) -> float | None:
        entry = pairs.get(str(row))
        if not isinstance(entry, dict) or str(column) not in entry:
            return None
        return float(entry[str(column)])

    values: list[float] = []
    for index, name in enumerate(order):
        if name not in binder:
            continue
        both = [v for v in (cell(index, peptide), cell(peptide, index)) if v is not None]
        if both:
            values.append(sum(both) / len(both))
    return sum(values) / len(values) if values else None


def read_confidence(run: Path) -> dict[str, object] | None:
    found = sorted(run.glob("confidence_*_model_0.json")) or sorted(
        run.glob("confidence_*.json")
    )
    return json.loads(found[0].read_text()) if found else None


@dataclass(frozen=True)
class Screen:
    """Every prediction read, plus the denominators needed to report what was not read."""

    runs: dict[str, dict[str, dict[str, object]]]
    directories: int
    unnamed: int
    unscored: int


def collect(predictions: Path) -> Screen:
    runs: dict[str, dict[str, dict[str, object]]] = {}
    directories = unnamed = unscored = 0
    for run in sorted(p for p in predictions.iterdir() if p.is_dir()):
        directories += 1
        parsed = parse_run_name(run.name)
        if parsed is None:
            unnamed += 1
            continue
        design, condition = parsed
        confidence = read_confidence(run)
        if confidence is None:
            unscored += 1
            continue
        runs.setdefault(design, {})[condition] = confidence
    return Screen(runs=runs, directories=directories, unnamed=unnamed, unscored=unscored)


def dependence_of(
    design: str,
    conditions: dict[str, dict[str, object]],
    chains: dict[str, dict[str, str]],
) -> Dependence | None:
    if COMPLEX not in conditions or APO not in conditions:
        return None
    complexed = conditions[COMPLEX]
    wildtype = conditions.get(WILDTYPE)
    order = chain_ids(design, tuple(chains[COMPLEX]))
    return Dependence(
        design=design,
        complex_iptm=float(complexed["iptm"]),
        apo_iptm=float(conditions[APO]["iptm"]),
        wildtype_iptm=None if wildtype is None else float(wildtype["iptm"]),
        complex_peptide_pair=peptide_pair(complexed, order, binder_ids(design)),
        complex_plddt=float(complexed["complex_plddt"]),
    )


@app.command()
def main(
    out_dir: Path = typer.Option(Path("cofold_out"), help="boltz --out_dir"),
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    out: Path = typer.Option(Path("cofold_dependence.csv")),
    top: int = typer.Option(25, help="how many to print"),
) -> None:
    predictions = next(out_dir.glob("boltz_results_*/predictions"), None)
    if predictions is None:
        raise typer.BadParameter(f"no boltz_results_*/predictions under {out_dir}")

    screen = collect(predictions)
    chains = conditions_from(targets)
    rows = [
        d
        for d in (dependence_of(k, v, chains) for k, v in screen.runs.items())
        if d is not None
    ]
    if not rows:
        raise typer.BadParameter(f"no design has both {COMPLEX} and {APO}")

    control = next((r for r in rows if r.design == CONTROL), None)
    designs = [r for r in rows if r.design != CONTROL]

    if control is None:
        typer.echo("NO CONTROL: every number below is uncalibrated and must be reported as such")
    else:
        wt = (
            f", wild type {control.wildtype_iptm:.3f}"
            if control.wildtype_iptm is not None
            else ""
        )
        typer.echo(
            f"control 7OW6 TCR: complex {control.complex_iptm:.3f}"
            f" apo {control.apo_iptm:.3f} -> drop {control.drop:+.3f}{wt}"
        )

    # Every denominator, because a design whose apo run failed drops out of the ranking and
    # would otherwise leave no trace in the one script that decides which 20 are submitted.
    scored = sum(len(v) for v in screen.runs.values())
    partial = sorted(k for k, v in screen.runs.items() if not {COMPLEX, APO} <= set(v))
    typer.echo(
        f"\n{screen.directories} prediction directories:"
        f" {scored} scored, {screen.unscored} with no confidence file,"
        f" {screen.unnamed} unrecognised names"
    )
    typer.echo(f"{len(rows)} designs have both conditions ({len(designs)} designs + control)")
    if partial:
        typer.echo(f"INCOMPLETE, excluded from ranking: {', '.join(partial)}")

    typer.echo(f"\n  {'design':20s} {'g12d':>6s} {'apo':>6s} {'drop':>7s} {'pair':>6s} {'plddt':>6s}")
    for row in sorted(designs, key=lambda r: -r.drop)[:top]:
        pair = "    --" if row.complex_peptide_pair is None else f"{row.complex_peptide_pair:6.3f}"
        typer.echo(
            f"  {row.design:20s} {row.complex_iptm:6.3f} {row.apo_iptm:6.3f}"
            f" {row.drop:+7.3f} {pair} {row.complex_plddt:6.3f}"
        )

    drops = sorted(r.drop for r in designs)
    typer.echo(
        f"\ndrop: min={drops[0]:+.3f} median={drops[len(drops) // 2]:+.3f} max={drops[-1]:+.3f}"
    )
    if control is not None:
        clearing = [r for r in designs if r.drop >= control.drop]
        typer.echo(
            f"{len(clearing)} of {len(designs)} lose at least as much confidence as the"
            f" receptor does ({control.drop:+.3f}) when the peptide is removed"
        )
    typer.echo(wrote(out, write_rows(rows, out)))


if __name__ == "__main__":
    app()
