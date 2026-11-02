"""Select the designs worth co-folding, and write the Boltz specs for them.

The selection lives here rather than in a notebook so the funnel is reproducible: the gates'
own CSVs go in, a manifest recording every design's fate comes out, and the specs written are
exactly the ones the manifest says survived.

Three conditions were prepared, but they do not all earn the same compute. The G12D-versus-wild-
type comparison is the first ranking tier, and this campaign's own measurement says it is
unreachable: the two complexes differ by a 0.5-0.85 A backbone shift, below the coordinate error
of Boltz-2 and Chai-1. Spending a prediction per design on a comparison that cannot resolve
would buy a number nobody may rely on. Peptide-dependence can be resolved -- the peptide-free
receptor is a different molecule, not a sub-angstrom rearrangement of one -- so the screen runs
`complex_g12d` against `apo_hla` on every survivor. `complex_wt` is written for the control only;
no design gets one from this screen, and a finalist's wild-type run comes from a later
invocation with `--conditions`, where it is reported rather than gated.

The 7OW6 receptor is written for all three conditions regardless. Without it an interface
confidence on this epitope has no scale: it is the only molecule known to do the job.

This module also owns the naming vocabulary -- run names, condition names and chain ids -- so
that `cofold_scores` imports it instead of recovering it from the YAML written here. The chain
order in a spec is ours by construction, and parsing generated output to rediscover it is a
second representation of something one import settles.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import typer
from common import rows_by, write_rows, wrote
from fold_novelty import FOLD_VERDICTS
from novelty import read_designs

app = typer.Typer(add_completion=False)

# The 7OW6 receptor's own peptide share, measured by engagement.py --control. A threshold, not
# an aspiration: the peptide is 1.2% of the complex's accessible surface and the one binder
# known to discriminate these peptides puts 55 of its 275 contacts on it. Only comparable to a
# share measured at the same contact cutoff, which engagement.py records per row.
CONTROL_PEPTIDE_SHARE = 55 / (55 + 220)

FOLD_REJECTS = ("copy_of_known", "ig_sandwich")
# Severity derived from the producer's ordering rather than restated here.
FOLD_SEVERITY = {verdict: len(FOLD_VERDICTS) - i for i, verdict in enumerate(FOLD_VERDICTS)}

COMPLEX = "complex_g12d"
APO = "apo_hla"
WILDTYPE = "complex_wt"
CONTROL = "control_tcr"

SCREEN_CONDITIONS = (COMPLEX, APO)

SEPARATOR = "__"
BINDER_CHAIN = "X"
CONTROL_BINDER_CHAINS = ("D", "E")
PEPTIDE_CHAIN = "C"


def run_name(design: str, condition: str) -> str:
    return f"{design}{SEPARATOR}{condition}"


def parse_run_name(name: str) -> tuple[str, str] | None:
    """Split a prediction directory name back into design and condition.

    Partitioned from the right because a condition name contains no separator while a design
    name might; tnf-binder recorded a funnel silently shrinking when the same split was taken
    from the left.
    """
    design, separator, condition = name.rpartition(SEPARATOR)
    return (design, condition) if separator else None


def binder_ids(design: str) -> tuple[str, ...]:
    return CONTROL_BINDER_CHAINS if design == CONTROL else (BINDER_CHAIN,)


def chain_ids(design: str, condition_chains: tuple[str, ...]) -> tuple[str, ...]:
    """Chain ids in the order the spec lists them, which is how Boltz indexes its pair matrix."""
    return binder_ids(design) + tuple(condition_chains)


def predictions_dir(out_dir: Path) -> Path:
    """Boltz's prediction directory inside an --out_dir. The layout is ours to know once."""
    found = next(out_dir.glob("boltz_results_*/predictions"), None)
    if found is None:
        raise typer.BadParameter(f"no boltz_results_*/predictions under {out_dir}")
    return found


def conditions_from(targets: Path) -> dict[str, dict[str, str]]:
    """Each co-fold condition's chains, in spec order, from prepare_target.py's target file."""
    return {entry["name"]: entry["chains"] for entry in json.loads(targets.read_text())}


@dataclass(frozen=True)
class Fate:
    design: str
    variant: str
    copied_fraction: float
    fold_verdict: str
    peptide_share: float
    hotspots_touched: int
    hotspots_nominated: int
    mutation_contacts: int
    verdict: str

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "variant": self.variant,
            "copied_fraction": self.copied_fraction,
            "fold_verdict": self.fold_verdict,
            "peptide_share": self.peptide_share,
            "hotspots_touched": self.hotspots_touched,
            "hotspots_nominated": self.hotspots_nominated,
            "mutation_contacts": self.mutation_contacts,
            "verdict": self.verdict,
        }


def worst_fold(fold_csv: Path) -> dict[str, str]:
    """Each design's most severe fold verdict across the databases searched."""
    worst: dict[str, str] = {}
    with fold_csv.open() as handle:
        for row in csv.DictReader(handle):
            design, verdict = row["design"], row["verdict"]
            if design not in worst or FOLD_SEVERITY[verdict] > FOLD_SEVERITY[worst[design]]:
                worst[design] = verdict
    return worst


def judge(design: str, nov: dict[str, str], fold: str, eng: dict[str, str]) -> Fate:
    share = float(eng["peptide_share"])
    touched, nominated = int(eng["hotspots_touched"]), int(eng["hotspots_nominated"])
    if nov["novel"] != "True":
        verdict = "rejected_sequence"
    elif fold in FOLD_REJECTS:
        verdict = f"rejected_{fold}"
    elif share < CONTROL_PEPTIDE_SHARE:
        verdict = "rejected_groove_gripper"
    elif touched < nominated:
        verdict = "rejected_missed_hotspot"
    else:
        verdict = "selected"
    return Fate(
        design=design,
        variant=eng["variant"],
        copied_fraction=float(nov["copied_fraction"]),
        fold_verdict=fold,
        peptide_share=share,
        hotspots_touched=touched,
        hotspots_nominated=nominated,
        mutation_contacts=int(eng["mutation_contacts"]),
        verdict=verdict,
    )


def design_sequences(gen_out: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    for metrics in sorted(gen_out.glob("*/final_ranked_designs/all_designs_metrics.csv")):
        found.update(read_designs(metrics, "id", "designed_sequence"))
    if not found:
        raise typer.BadParameter(f"no all_designs_metrics.csv under {gen_out}")
    return found


def spec_text(chains: list[tuple[str, str]]) -> str:
    """A Boltz 2 co-folding spec. `msa: empty` because --use_msa_server is a flag, not a value,
    and a design has no homologues to align anyway."""
    blocks = [
        f"  - protein:\n      id: {chain_id}\n      sequence: {sequence}\n      msa: empty"
        for chain_id, sequence in chains
    ]
    return "version: 1\nsequences:\n" + "\n".join(blocks) + "\n"


def write_spec(out: Path, name: str, chains: list[tuple[str, str]]) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.yaml"
    path.write_text(spec_text(chains))
    return path


Entry = tuple[str, list[tuple[str, str]]]


def design_entries(
    designs: Sequence[str],
    sequences: Mapping[str, str],
    conditions: Mapping[str, Mapping[str, str]],
    condition_names: Sequence[str],
) -> list[Entry]:
    """One spec per design per condition: the binder as chain X, then the condition's chains."""
    return [
        (
            run_name(design, condition),
            [(BINDER_CHAIN, sequences[design]), *conditions[condition].items()],
        )
        for design in designs
        for condition in condition_names
    ]


def control_entries(
    tcr: Mapping[str, str],
    conditions: Mapping[str, Mapping[str, str]],
    condition_names: Sequence[str] | None = None,
) -> list[Entry]:
    """The 7OW6 receptor through the given conditions, so every pass carries its own scale.

    Defaults to every condition: the screen runs the control through all three even though the
    designs only get two, because the wild-type number is only interpretable next to it.
    """
    return [
        (run_name(CONTROL, condition), [*tcr.items(), *conditions[condition].items()])
        for condition in (tuple(conditions) if condition_names is None else condition_names)
    ]


def write_entries(out: Path, entries: Sequence[Entry]) -> int:
    for name, chains in entries:
        write_spec(out, name, chains)
    return len(entries)


@app.command()
def main(
    novelty: Path = typer.Option(Path("."), help="directory holding novelty_<variant>.csv"),
    fold: Path = typer.Option(Path("fold_novelty_cohort.csv")),
    engagement: Path = typer.Option(Path("engagement_cohort.csv")),
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    controls: Path = typer.Option(Path("targets/controls.json")),
    gen_out: Path = typer.Option(Path("gen_out")),
    out: Path = typer.Option(Path("cofold_specs")),
    manifest: Path = typer.Option(Path("cofold_manifest.csv")),
) -> None:
    nov = {
        design: row
        for path in sorted(novelty.glob("novelty_*.csv"))
        for design, row in rows_by(path, "design").items()
    }
    if not nov:
        raise typer.BadParameter(f"no novelty_*.csv in {novelty}")
    fold_verdicts = worst_fold(fold)
    eng = rows_by(engagement, "design")
    sequences = design_sequences(gen_out)

    fates = [
        judge(design, row, fold_verdicts.get(design, "distinct_fold"), eng[design])
        for design, row in nov.items()
        if design in eng
    ]
    selected = [f for f in fates if f.verdict == "selected"]

    conditions = conditions_from(targets)
    if missing := [name for name in SCREEN_CONDITIONS if name not in conditions]:
        raise typer.BadParameter(f"{targets} has no condition(s) {missing}")
    tcr = json.loads(controls.read_text())["tcr"]

    # Resolve every spec before writing any of it: a selected design with no sequence in
    # gen_out used to raise part-way through the loop, leaving a half-written spec directory
    # and -- because the manifest was written last -- no record of the funnel at all.
    if absent := sorted(f.design for f in selected if f.design not in sequences):
        raise typer.BadParameter(f"{gen_out} has no sequence for selected design(s) {absent}")

    specs = design_entries(
        [f.design for f in selected], sequences, conditions, SCREEN_CONDITIONS
    ) + control_entries(tcr, conditions)

    typer.echo(wrote(manifest, write_rows(fates, manifest)))
    written = write_entries(out, specs)

    tally = Counter(f.verdict for f in fates)
    typer.echo("funnel over every generated design:")
    for verdict, count in tally.most_common():
        typer.echo(f"  {verdict:28s} {count:4d}")
    typer.echo(
        f"selected {len(selected)} of {len(fates)};"
        f" by variant {dict(Counter(f.variant for f in selected))}"
    )
    typer.echo(f"wrote {written} specs to {out}")


if __name__ == "__main__":
    app()
