"""The co-fold passes that run after the screen has ranked the survivors.

**Repeated sampling.** The screen drew one diffusion sample per condition, so a design's drop
was a difference of two single draws and carried the variance of both. `select-noise` wrote
specs for a subset spanning that ranking to be re-run at five samples, and the probe answered
its question in the worst way available: two of the pilot's top four reversed sign. The whole
cohort is therefore re-run at five samples, `report-spread` reads it, and the one-sample screen
is retrospectively a pilot. `select-noise` and `spanning` are kept as the record of how that
was established; nothing downstream reads their output any more.

**One significance rule, defined here.** `report-spread` and `submission.py` both have to say
whether a difference is real, and they must say it the same way or the table a human reads
disagrees with the file that gets uploaded. `SIGMA` and `comparison` below are that rule, and
they are here rather than in the submission stage because this is the only module that sees the
samples -- in particular it is the only place that knows `n`, and `n` is what makes the rule
correct. `drop_sd` is the spread of *one* difference; the quantity being tested is a mean of
`n`, whose standard error is `drop_sem = drop_sd / sqrt(n)`. Testing a mean against the
per-draw spread is a sqrt(n)-times-too-strict gate, which under-reports real drops.

This is not a fresh discovery: `touchstone.boltz` already records the measurement that best and
mean differ by up to 0.35 within a single design, and that ranking on best-of-n reorders the
field. What is new here is that the ranked quantity is a *difference* of two sampled scores, so
it carries both conditions' variance.

**Wild type.** `select-wildtype` writes the `complex_wt` condition for the finalists only. It is
reported, never gated, and `report-replicates` is why: paired across eight independent
invocations the control scores wild type *above* G12D by +0.0027 +- 0.0006, in seven of the
eight. The one receptor known to discriminate these peptides is scored here as preferring the
peptide it does not recognise -- so on the only molecule that can calibrate ranking tier 1,
this predictor's sign is not merely unresolved but inverted, and no selectivity a design shows
in this column could mean anything. The column exists so the claim is auditable rather than
asserted.

**The bar is reproducible, and the within-run error bar is already honest.** Two 5-sample runs
of a byte-identical control spec gave drops of +0.104 (26-spec probe) and +0.150 (159-spec
cohort), so the bar the campaign is calibrated against appeared to move by four times its own
standard error. More samples inside one invocation cannot tell a between-run component from a
bad draw -- Boltz seeds once per process, so a run's samples share that state -- which is why
`report-replicates` reads *separate* invocations. Eight of them put the drop at +0.101 with a
spread of 0.004, against the 0.007 each run reports internally; a permutation test over the
pooled samples puts that spread at the 4th percentile of the no-run-effect null. So there is no
resolved between-run component, `drop_sem` needs no widening, and the cohort's +0.150 is not
run-to-run noise but something specific to that one invocation -- which is why the cohort pass
is replicated in full rather than this bar being swapped in for it. `run_effect_sd` is the
decomposition that says so, and it subtracts: `Var(run means)` already contains `within_sem**2`.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, stdev

import typer
from cofold_scores import collect, sample_confidences
from cofold_specs import (
    APO,
    COMPLEX,
    CONTROL,
    SCREEN_CONDITIONS,
    WILDTYPE,
    conditions_from,
    control_entries,
    design_entries,
    design_sequences,
    find_predictions,
    predictions_dir,
    write_entries,
)
from common import rounded, rows_by, write_rows, wrote

app = typer.Typer(add_completion=False)

# How many standard errors two means must be apart before the difference is called a difference.
# Two, so the claim is roughly 95% rather than "the means differ", which at these spreads they
# always do.
SIGMA = 2.0


def comparison(mean_a: float, sem_a: float, mean_b: float, sem_b: float) -> str:
    """Whether mean_a sits above, tied with, or below mean_b, given each one's standard error.

    The one place the campaign decides what counts as a difference. Pass `0.0, 0.0` for the
    second pair to test a single mean against zero.
    """
    margin = SIGMA * math.hypot(sem_a, sem_b)
    difference = mean_a - mean_b
    if difference > margin:
        return "above"
    return "below" if difference < -margin else "tied"


@dataclass(frozen=True)
class Paired:
    """A quantity measured once per independent Boltz invocation of the same spec.

    `Spread` sees only the samples inside one invocation, and Boltz seeds once per process, so
    a run's samples share that state: its standard error is the diffusion's spread around
    whatever that process drew, not the spread of the drawn value across processes. Only
    separate invocations reach the second component. `Paired` is every statistic that follows
    from them, used for the control's bar and for the wild-type difference alike, because both
    are one number per run.
    """

    label: str
    values: tuple[float, ...]

    @property
    def runs(self) -> int:
        return len(self.values)

    @property
    def mean(self) -> float:
        return mean(self.values)

    @property
    def between_sd(self) -> float:
        """Spread of one invocation's value around the pooled mean."""
        return stdev(self.values)

    @property
    def sem(self) -> float:
        """Standard error of the pooled mean, which averages the between-run spread down."""
        return self.between_sd / math.sqrt(self.runs)

    @property
    def verdict(self) -> str:
        return comparison(self.mean, self.sem, 0.0, 0.0)

    @property
    def row(self) -> dict[str, object]:
        return {
            "quantity": self.label,
            "runs": self.runs,
            "mean": round(self.mean, 4),
            "between_sd": round(self.between_sd, 4),
            "sem": round(self.sem, 4),
            "verdict": self.verdict,
        }


def run_effect_sd(between_sd: float, within_sem: float) -> float:
    """The part of the run-to-run spread that sampling inside a run does not already explain.

    A run's mean is its own draw plus that run's sampling error, so `Var(run means)` already
    *contains* `within_sem**2`: the two subtract, they do not add in quadrature. Measured on the
    control, `between_sd` 0.004 against `within_sem` 0.007 leaves a negative variance, and a
    permutation test over the pooled samples puts the observed spread at the 4th percentile of
    the no-run-effect null -- so the drop has no resolved between-run component, and one
    invocation's own standard error is already the honest error bar (if anything conservative).

    Clamped at zero rather than hidden, because "not resolved" is the finding.
    """
    return math.sqrt(max(0.0, between_sd**2 - within_sem**2))


def spanning(designs: list[tuple[str, float]], count: int, head: int = 0) -> list[str]:
    """The top `head` of the drop ranking, plus the rest spread evenly over it.

    Spread over the *ranking* rather than the value, so the subset is not crowded into whatever
    part of the range happens to be densest. The head exists because the contenders are where
    the noise estimate has to settle an actual decision: whether the designs at the top of the
    list are really ahead of the ones below them, or tied.
    """
    if count >= len(designs):
        return [name for name, _ in designs]
    # One slot always goes to the far end, or the subset would span nothing -- which is what a
    # head at least as large as the count used to produce, silently.
    head = min(head, count - 1)
    last, remaining = len(designs) - 1, count - head
    picked = set(range(head))
    if remaining == 1:
        # Only reachable when head == count - 1. With a head to span away from, the one free
        # slot is the far end; with no head at all, count is 1 and the top is what to measure.
        picked.add(last if head else 0)
    else:
        step = (last - head) / (remaining - 1)
        picked |= {round(head + i * step) for i in range(remaining)}
    return [designs[i][0] for i in sorted(picked)]


def ranked_designs(dependence: Path) -> list[tuple[str, float]]:
    """Designs by peptide-dependence drop, best first, with the control left out."""
    ranked = sorted(
        (
            (design, float(row["peptide_dependence_drop"]))
            for design, row in rows_by(dependence, "design").items()
            if design != CONTROL
        ),
        key=lambda pair: -pair[1],
    )
    if not ranked:
        raise typer.BadParameter(f"{dependence} holds no designs")
    return ranked


@app.command()
def select_noise(
    dependence: Path = typer.Option(Path("cofold_dependence.csv")),
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    gen_out: Path = typer.Option(Path("gen_out")),
    controls: Path = typer.Option(Path("targets/controls.json")),
    out: Path = typer.Option(Path("cofold_noise_specs")),
    count: int = typer.Option(12, help="designs to re-run, spanning the drop ranking"),
    head: int = typer.Option(4, help="how many of the top-ranked to include regardless"),
) -> None:
    ranked = ranked_designs(dependence)
    chosen = spanning(ranked, count, head)
    drops = dict(ranked)
    conditions = conditions_from(targets)
    sequences = design_sequences(gen_out)
    tcr = json.loads(controls.read_text())["tcr"]

    written = write_entries(
        out,
        design_entries(chosen, sequences, conditions, SCREEN_CONDITIONS)
        + control_entries(tcr, conditions, SCREEN_CONDITIONS),
    )

    typer.echo(f"spanning subset of {len(chosen)} designs plus the control:")
    for design in chosen:
        typer.echo(f"  {design:20s} screen drop {drops[design]:+.3f}")
    typer.echo(f"wrote {written} specs to {out}")


@dataclass(frozen=True)
class Spread:
    """One design's repeated ipTM draws, and the statistics that follow from them."""

    design: str
    complexed: tuple[float, ...]
    apo: tuple[float, ...]

    @property
    def samples(self) -> str:
        """Both denominators: the two conditions are sampled independently."""
        return f"{len(self.complexed)}/{len(self.apo)}"

    @property
    def complex_mean(self) -> float:
        return mean(self.complexed)

    @property
    def complex_sd(self) -> float:
        return stdev(self.complexed)

    @property
    def apo_mean(self) -> float:
        return mean(self.apo)

    @property
    def apo_sd(self) -> float:
        return stdev(self.apo)

    @property
    def drop_mean(self) -> float:
        return self.complex_mean - self.apo_mean

    @property
    def drop_sd(self) -> float:
        """Spread of a difference of two independent draws, one per condition."""
        return math.hypot(self.complex_sd, self.apo_sd)

    @property
    def drop_sem(self) -> float:
        """Standard error of `drop_mean` -- the quantity actually compared.

        Each condition's mean carries `sd / sqrt(n)`, and the drop is their difference, so the
        two combine in quadrature. This, not `drop_sd`, is the denominator of the tie rule: the
        conditions are sampled independently and may carry different `n`, which is why the
        division happens per condition rather than once at the end.
        """
        return math.hypot(
            self.complex_sd / math.sqrt(len(self.complexed)),
            self.apo_sd / math.sqrt(len(self.apo)),
        )

    @property
    def drop_significant(self) -> bool:
        """Whether this drop clears its own sampling noise."""
        return comparison(self.drop_mean, self.drop_sem, 0.0, 0.0) == "above"

    @property
    def drop_min(self) -> float:
        """Widest the drop gets across every pairing of the two conditions' samples.

        Not a paired statistic: no complex draw belongs to any particular apo draw, so the
        extremes come from the best complex against the worst apo and vice versa. The range
        widens with the number of samples, so it is audit material rather than a quantity to
        compare between designs sampled differently -- drop_sd is what the tie rule uses.
        """
        return min(self.complexed) - max(self.apo)

    @property
    def drop_max(self) -> float:
        return max(self.complexed) - min(self.apo)

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "samples": self.samples,
            "iptm_complex_mean": round(self.complex_mean, 4),
            "iptm_complex_sd": round(self.complex_sd, 4),
            "iptm_apo_mean": round(self.apo_mean, 4),
            "iptm_apo_sd": round(self.apo_sd, 4),
            "drop_mean": round(self.drop_mean, 4),
            "drop_sd": round(self.drop_sd, 4),
            "drop_sem": round(self.drop_sem, 4),
            "drop_significant": self.drop_significant,
            "drop_min": round(self.drop_min, 4),
            "drop_max": round(self.drop_max, 4),
        }


def sample_iptms(run: Path) -> tuple[float, ...] | None:
    """Every diffusion sample's ipTM for one prediction, or None if the run produced none."""
    values = tuple(float(c["iptm"]) for c in sample_confidences(run))
    return values or None


def spread_of(design: str, conditions: dict[str, tuple[float, ...]]) -> Spread | None:
    complexed, apo = conditions.get(COMPLEX, ()), conditions.get(APO, ())
    if len(complexed) < 2 or len(apo) < 2:
        return None
    return Spread(design=design, complexed=complexed, apo=apo)


@app.command()
def select_control(
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    controls: Path = typer.Option(Path("targets/controls.json")),
    out: Path = typer.Option(Path("cofold_control_specs")),
) -> None:
    """The 7OW6 receptor alone, through every condition, for the replicated bar measurement.

    Its own spec directory rather than a `cp` out of `cofold_specs`, so the one pass that
    calibrates the campaign reads specs written by the module that owns the naming.
    """
    tcr = json.loads(controls.read_text())["tcr"]
    written = write_entries(out, control_entries(tcr, conditions_from(targets)))
    typer.echo(f"wrote {written} control specs to {out}")


@app.command()
def select_wildtype(
    dependence: Path = typer.Option(Path("cofold_dependence.csv")),
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    gen_out: Path = typer.Option(Path("gen_out")),
    out: Path = typer.Option(Path("cofold_wt_specs")),
    count: int = typer.Option(20, help="finalists, which is the submission cap"),
) -> None:
    chosen = ranked_designs(dependence)[:count]
    finalists = [name for name, _ in chosen]
    conditions = conditions_from(targets)
    if WILDTYPE not in conditions:
        raise typer.BadParameter(f"{targets} has no condition {WILDTYPE}")

    written = write_entries(
        out,
        design_entries(finalists, design_sequences(gen_out), conditions, (WILDTYPE,)),
    )
    typer.echo(
        f"{len(chosen)} finalists by drop,"
        f" {chosen[0][1]:+.3f} down to {chosen[-1][1]:+.3f}"
    )
    typer.echo(f"wrote {written} specs to {out}")


def invocation_predictions(out_dir: Path) -> list[Path]:
    """Each independent invocation's predictions directory under `out_dir`.

    A replicated pass writes one sub-directory per invocation; a single pass writes predictions
    straight into `out_dir`. Both are the same thing here -- a list of separate processes -- so
    either shape is accepted rather than making the caller know which it was handed. A
    sub-directory holding no predictions is not an invocation and is skipped, so a stray
    sibling cannot inflate the denominator.
    """
    if (found := find_predictions(out_dir)) is not None:
        return [found]
    nested = (find_predictions(path) for path in sorted(out_dir.iterdir()) if path.is_dir())
    if not (invocations := [path for path in nested if path is not None]):
        raise typer.BadParameter(f"{out_dir} holds no predictions and no replicate directories")
    return invocations


@dataclass(frozen=True)
class RunRow:
    """One design in one invocation: its within-run statistics, labelled by which run."""

    replicate: str
    spread: Spread
    wildtype_mean: float | None

    @property
    def row(self) -> dict[str, object]:
        return {"replicate": self.replicate, **self.spread.row, "iptm_wt_mean": rounded(self.wildtype_mean)}


@dataclass(frozen=True)
class Walk:
    """Every design-run read, plus the denominators needed to report what was not read."""

    rows: tuple[RunRow, ...]
    directories: int
    unscored: int
    unnamed: int


def run_rows(out_dirs: Sequence[Path]) -> Walk:
    """Every design's within-run statistics from every invocation, in invocation order.

    `collect` counts what it could not score as well as what it could, and those counts are
    carried out of here rather than dropped: a replicate missing two of its predictions must
    not read like one that produced them all.
    """
    rows: list[RunRow] = []
    directories = unscored = unnamed = 0
    for out_dir in out_dirs:
        for predictions in invocation_predictions(out_dir):
            screen = collect(predictions, sample_iptms)
            directories += screen.directories
            unscored += screen.unscored
            unnamed += screen.unnamed
            for design, conditions in sorted(screen.runs.items()):
                if (spread := spread_of(design, conditions)) is None:
                    continue
                wildtype = conditions.get(WILDTYPE)
                rows.append(
                    RunRow(
                        replicate=predictions.parent.parent.name,
                        spread=spread,
                        wildtype_mean=mean(wildtype) if wildtype else None,
                    )
                )
    if not rows:
        raise typer.BadParameter(f"no design in {[str(d) for d in out_dirs]} has repeated samples")
    return Walk(tuple(rows), directories=directories, unscored=unscored, unnamed=unnamed)


def repeated(rows: Sequence[RunRow]) -> dict[str, list[RunRow]]:
    """Designs measured by more than one invocation -- the only ones this pass can speak for."""
    by_design: dict[str, list[RunRow]] = {}
    for row in rows:
        by_design.setdefault(row.spread.design, []).append(row)
    return {design: runs for design, runs in by_design.items() if len(runs) > 1}


def pooled_between_sd(sds: Sequence[float]) -> float:
    """One between-run SD for the whole cohort, pooling each design's own across designs.

    Variances pool, not standard deviations. At two invocations per design an individual
    estimate carries about 70% relative error, so a per-design number would be noise dressed as
    a property of the design; pooling buys roughly one degree of freedom per design.
    """
    return math.sqrt(mean(sd**2 for sd in sds))


def paired_drops(by_design: Mapping[str, list[RunRow]]) -> dict[str, Paired]:
    """Each design's drop across invocations, which is what every statistic below is built on."""
    return {
        design: Paired(design, tuple(row.spread.drop_mean for row in runs))
        for design, runs in by_design.items()
    }


def echo_replicate_table(by_design: Mapping[str, list[RunRow]], paired: Mapping[str, Paired]) -> None:
    typer.echo(
        f"\n  {'design':20s} {'runs':>4s} {'drop':>8s} {'between':>8s}"
        f" {'within':>7s} {'effect':>7s} {'vs 0':>6s}  per-run drops"
    )
    for design, drops in sorted(paired.items(), key=lambda kv: -kv[1].mean):
        within = mean(row.spread.drop_sem for row in by_design[design])
        typer.echo(
            f"  {design:20s} {drops.runs:4d} {drops.mean:+8.3f} {drops.between_sd:8.3f}"
            f" {within:7.3f} {run_effect_sd(drops.between_sd, within):7.3f} {drops.verdict:>6s}"
            f"  {' '.join(f'{v:+.3f}' for v in drops.values)}"
        )


@app.command()
def report_replicates(
    out_dir: list[Path] = typer.Option(
        [Path("cofold_ctrl_out")],
        help="repeatable: a directory of replicate sub-directories, or one pass's --out_dir",
    ),
    out: Path = typer.Option(Path("cofold_replicates.csv")),
) -> None:
    walk = run_rows(out_dir)
    by_design = repeated(walk.rows)
    if not by_design:
        raise typer.BadParameter(
            f"no design appears in more than one of {[str(d) for d in out_dir]};"
            " a between-run spread needs at least two independent invocations"
        )
    invocations = len({row.replicate for row in walk.rows})
    typer.echo(
        f"{walk.directories} prediction directories over {invocations} invocations:"
        f" {len(walk.rows)} design-runs scored, {walk.unscored} with no samples,"
        f" {walk.unnamed} unrecognised names"
    )
    singles = len({row.spread.design for row in walk.rows}) - len(by_design)
    typer.echo(
        f"{len(by_design)} designs measured more than once;"
        f" {singles} appear in one invocation only and are excluded"
    )

    paired = paired_drops(by_design)
    echo_replicate_table(by_design, paired)
    echo_bar(by_design, paired)
    echo_wildtype(walk.rows)
    typer.echo(wrote(out, write_rows(walk.rows, out)))


def echo_bar(by_design: Mapping[str, list[RunRow]], paired: Mapping[str, Paired]) -> None:
    """The control's drop, and whether the spread across invocations needs its own error term."""
    if (control := paired.get(CONTROL)) is None:
        typer.echo(
            f"\nNO CONTROL: {CONTROL} is not repeated here, so this pass calibrates nothing."
            " Nothing downstream may rank against it; re-run with the control's specs included."
        )
        return
    within = mean(row.spread.drop_sem for row in by_design[CONTROL])
    effect = run_effect_sd(control.between_sd, within)
    typer.echo(
        f"\nBAR: the control's drop is {control.mean:+.3f} over {control.runs} invocations,"
        f" spread {control.between_sd:.3f} against a within-run error of {within:.3f}"
    )
    typer.echo(
        f"run effect {effect:.3f}:"
        + (
            " not resolved above within-run sampling, so one invocation's own standard error is"
            " already the honest bar and needs no between-run widening"
            if effect == 0.0
            else " a design measured once carries this on top of its own standard error"
        )
    )

    designs = [p for design, p in paired.items() if design != CONTROL]
    if designs:
        typer.echo(
            f"pooled over {len(designs)} designs: between-run SD"
            f" {pooled_between_sd([p.between_sd for p in designs]):.3f} against a within-run"
            f" error of {mean(row.spread.drop_sem for runs in by_design.values() for row in runs):.3f}"
        )


def echo_wildtype(rows: Sequence[RunRow]) -> None:
    """Ranking tier 1, paired within each invocation and measured on the control alone.

    Control-only by construction rather than by whichever specs the pass happened to hold: the
    claim is about the one receptor known to discriminate these peptides, and pooling designs
    into it would change what the number means while printing the same sentence.
    """
    paired = [
        (row.wildtype_mean, row.spread.complex_mean)
        for row in rows
        if row.spread.design == CONTROL and row.wildtype_mean is not None
    ]
    if len(paired) < 2:
        return
    effect = Paired(f"{WILDTYPE}_minus_{COMPLEX}", tuple(wt - cx for wt, cx in paired))
    typer.echo(
        f"\nTIER 1: wild type minus G12D is {effect.mean:+.4f} +- {effect.sem:.4f} over"
        f" {effect.runs} invocations -- {effect.verdict} zero."
        " Positive means this predictor scores the peptide the control does NOT recognise"
        " higher, so its sign carries no information about a design's selectivity"
    )


@app.command()
def report_spread(
    out_dir: Path = typer.Option(Path("cofold_full5_out"), help="boltz --out_dir"),
    out: Path = typer.Option(Path("cofold_spread.csv")),
) -> None:
    screen = collect(predictions_dir(out_dir), sample_iptms)
    typer.echo(
        f"{screen.directories} prediction directories:"
        f" {sum(len(v) for v in screen.runs.values())} scored,"
        f" {screen.unscored} with no samples, {screen.unnamed} unrecognised names"
    )

    rows = [s for s in (spread_of(k, v) for k, v in screen.runs.items()) if s is not None]
    if not rows:
        raise typer.BadParameter(f"no design in {out_dir} has repeated samples for both conditions")

    control = next((r for r in rows if r.design == CONTROL), None)
    designs = [r for r in rows if r.design != CONTROL]

    typer.echo(
        f"\n  {'design':20s} {'n':>5s} {'complex':>15s} {'apo':>15s}"
        f" {'drop':>8s} {'sem':>6s} {'>noise':>7s} {'vs ctrl':>8s} {'range':>16s}"
    )
    for row in sorted(rows, key=lambda r: -r.drop_mean):
        versus = "--" if control is None else comparison(
            row.drop_mean, row.drop_sem, control.drop_mean, control.drop_sem
        )
        typer.echo(
            f"  {row.design:20s} {row.samples:>5s}"
            f" {row.complex_mean:7.3f}+-{row.complex_sd:<6.3f}"
            f" {row.apo_mean:7.3f}+-{row.apo_sd:<6.3f}"
            f" {row.drop_mean:+8.3f} {row.drop_sem:6.3f}"
            f" {'yes' if row.drop_significant else 'no':>7s} {versus:>8s}"
            f" [{row.drop_min:+.3f},{row.drop_max:+.3f}]"
        )

    # The same rule the submission is built from, printed here so the table a human reads and
    # the file that gets uploaded cannot disagree about which designs are tied.
    significant = sum(r.drop_significant for r in designs)
    typer.echo(
        f"\n{significant} of {len(designs)} designs have a drop clearing its own standard error"
        f" at {SIGMA:g} sigma"
    )
    if control is None:
        typer.echo(
            f"NO CONTROL: {CONTROL} has no repeated samples here, so these drops have no scale."
            " Nothing downstream may rank them; re-run the pass with the control's specs."
        )
    else:
        verdicts = [
            comparison(r.drop_mean, r.drop_sem, control.drop_mean, control.drop_sem)
            for r in designs
        ]
        tally = {v: verdicts.count(v) for v in ("above", "tied", "below")}
        typer.echo(
            f"against the control's {control.drop_mean:+.3f} +- {control.drop_sem:.3f}"
            f" over {control.samples} samples: {tally['above']} above,"
            f" {tally['tied']} tied, {tally['below']} below -- its bar is a band, not a line"
        )
    typer.echo(wrote(out, write_rows(rows, out)))


if __name__ == "__main__":
    app()
