"""Pool BoltzGen's per-spec metrics into one table, labelled by epitope and band.

BoltzGen writes `<out>/<spec>/final_ranked_designs/all_designs_metrics.csv` with ~230
columns. Most are not useful here and a few are actively misleading, so this keeps a
declared subset and says why each is kept:

- `design_to_target_iptm` -- the generator's own interface confidence. Kept to be
  *compared against an independent model*, not to rank by: on 2,214 wet-lab-tested public
  designs ipTM's within-target AUROC for binding was 0.607, and the observed binding rate
  fell across most of its range. It is a screening floor, nothing more.
- `delta_sasa_refolded`, `plip_hbonds_refolded`, `plip_saltbridge_refolded` -- interface
  size and polar contacts, measured on the refolded complex.
- `HIS_fraction`, `ASP_fraction`, `GLU_fraction` -- this challenge is ranked on a pH switch
  first, so histidine content and the acidic context are first-class here.
- `GLU_fraction` doubles as the expression predictor that replicated at AUROC 0.724 on the
  public corpus, better than any structural score tested.
- `liability_UnpairedCys_count` -- 18 of 20 designs in the TNF submission carried an
  unpaired cysteine from one chassis, found only at packaging time.
- `num_filters_passed` / `pass_filters` -- BoltzGen's own verdict, recorded rather than
  trusted, so our triage and its triage can be compared.

Two failure modes this is built to refuse:

- **A renamed column reads as a blank one.** BoltzGen moves its confidence key names
  between versions -- `touchstone/src/touchstone/boltzgen.py` records the same thing and
  keeps them as fallback chains -- so every column here is a chain of candidate names and
  a spec that carries none of them is an error naming the column, not a silently empty
  field that the next stage ranks on.
- **"No designs" has three causes and one of them is a crash.** A missing metrics file, a
  metrics file with no rows, and a file whose every row has an empty sequence are
  different events. They are written to a `.specs.json` sidecar beside the pooled CSV so
  the funnel's denominators survive the terminal session that produced them.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import typer
from csvio import write_rows

# short name -> candidate BoltzGen column names, first one present wins.
KEEP: dict[str, tuple[str, ...]] = {
    "iptm": ("design_to_target_iptm", "ligand_iptm", "iptm"),
    "pae": ("min_design_to_target_pae",),
    "ptm": ("design_ptm", "ptm"),
    "dsasa": ("delta_sasa_refolded",),
    "hbonds": ("plip_hbonds_refolded",),
    "saltbridges": ("plip_saltbridge_refolded",),
    "his_frac": ("HIS_fraction",),
    "asp_frac": ("ASP_fraction",),
    "glu_frac": ("GLU_fraction",),
    "cys_frac": ("CYS_fraction",),
    "ala_frac": ("ALA_fraction",),
    "gly_frac": ("GLY_fraction",),
    "unpaired_cys": ("liability_UnpairedCys_count",),
    "liability": ("liability_score",),
    "filters_passed": ("num_filters_passed",),
    "boltzgen_pass": ("pass_filters",),
    "bb_rmsd": ("bb_rmsd",),
}

SEQUENCE_COLUMNS = ("designed_sequence", "designed_chain_sequence", "sequence")
ID_COLUMNS = ("id", "design_id", "name")

OUT_COLUMNS = ["name", "epitope", "band", "length", *KEEP, "sequence"]

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class SpecResult:
    """What one spec contributed, including the denominators when it contributed nothing."""

    spec: str
    state: str
    rows_read: int = 0
    dropped_no_sequence: int = 0
    designs: list[dict[str, object]] = field(default_factory=list)

    @property
    def note(self) -> str:
        if self.state == "ok":
            dropped = (f", {self.dropped_no_sequence} dropped with no sequence"
                       if self.dropped_no_sequence else "")
            return f"{len(self.designs):4d} designs{dropped}"
        return {
            "no_metrics_file": "   no metrics file - this spec was never scored",
            "no_rows": "   metrics file present but empty - the spec ran and produced nothing",
            "no_sequences": "   every row had an empty sequence - scored but unusable",
        }[self.state]

    def summary(self) -> dict[str, object]:
        return {
            "spec": self.spec,
            "state": self.state,
            "rows_read": self.rows_read,
            "dropped_no_sequence": self.dropped_no_sequence,
            "designs_kept": len(self.designs),
        }


def epitope_and_band(spec: str) -> tuple[str, str]:
    return (spec[:-2], "short") if spec.endswith("_s") else (spec, "long")


def resolve(fieldnames: list[str], wanted: dict[str, tuple[str, ...]]) -> dict[str, str]:
    """Map each short name to the column the file actually carries, or raise naming the gaps."""
    present = set(fieldnames)
    found = {short: next((c for c in candidates if c in present), None)
             for short, candidates in wanted.items()}
    if missing := {short: wanted[short] for short, column in found.items() if column is None}:
        raise typer.BadParameter(
            "the metrics file carries none of the candidate names for "
            + "; ".join(f"{short} ({'/'.join(names)})" for short, names in missing.items())
            + ". BoltzGen renames these between versions, so this is a version mismatch "
              "rather than a missing measurement -- adding the new name to KEEP is the fix, "
              "and ranking on a blank column is what this refuses to do."
        )
    return {short: column for short, column in found.items() if column is not None}


def rows_for_spec(directory: Path) -> SpecResult:
    metrics = directory / "final_ranked_designs" / "all_designs_metrics.csv"
    if not metrics.exists():
        return SpecResult(directory.name, "no_metrics_file")

    with metrics.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = resolve(list(reader.fieldnames or []), KEEP)
        identity = next((c for c in ID_COLUMNS if c in (reader.fieldnames or [])), None)
        sequences = next((c for c in SEQUENCE_COLUMNS if c in (reader.fieldnames or [])), None)
        rows = list(reader)

    if not rows:
        return SpecResult(directory.name, "no_rows")
    if sequences is None or identity is None:
        raise typer.BadParameter(
            f"{metrics}: no usable {'sequence' if sequences is None else 'id'} column "
            f"(tried {SEQUENCE_COLUMNS if sequences is None else ID_COLUMNS})"
        )

    epitope, band = epitope_and_band(directory.name)
    designs: list[dict[str, object]] = []
    for row in rows:
        sequence = (row.get(sequences) or "").strip().upper()
        if not sequence:
            continue
        designs.append({
            "name": f"{directory.name}_{row[identity]}",
            "epitope": epitope,
            "band": band,
            "length": len(sequence),
            "sequence": sequence,
            **{short: row.get(column, "") for short, column in columns.items()},
        })

    state = "ok" if designs else "no_sequences"
    return SpecResult(directory.name, state, len(rows), len(rows) - len(designs), designs)


@app.command()
def main(
    pilot_out: Annotated[Path, typer.Option(help="BoltzGen output root")],
    out: Annotated[Path, typer.Option(help="pooled per-design CSV")],
) -> None:
    directories = sorted(d for d in pilot_out.iterdir() if d.is_dir())
    if not directories:
        raise typer.BadParameter(f"no spec directories under {pilot_out}")

    results = []
    for directory in directories:
        result = rows_for_spec(directory)
        results.append(result)
        typer.echo(f"  {result.spec:16s} {result.note}")

    designs = [design for result in results for design in result.designs]
    if not designs:
        raise typer.BadParameter("no designs with sequences found in any spec")

    write_rows(out, designs, OUT_COLUMNS)
    sidecar = out.with_suffix(".specs.json")
    sidecar.write_text(json.dumps([r.summary() for r in results], indent=2) + "\n")

    typer.echo(f"\n{len(designs)} designs pooled from {len(results)} specs "
               f"({sum(1 for r in results if r.state == 'ok')} scored)")
    for epitope, count in sorted(Counter(str(d["epitope"]) for d in designs).items()):
        typer.echo(f"  {epitope:14s} {count:4d}")
    typer.echo(f"wrote {out} and {sidecar}")


if __name__ == "__main__":
    app()
