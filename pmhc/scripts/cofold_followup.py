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

**Which table the submission ranks on.** `report-spread` reads one pass; `pool-replicates`
reads several and writes the same column shape, so `submission.py` consumes either without
knowing which. That is deliberate: the measurement of record is whatever the best available
table is, and making the pooled one a *file* keeps it an artifact a reader can check rather
than arithmetic hidden inside the submission step. The columns mean different things between
the two tables -- see `Pooled` -- which is the price of the single format.

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
from common import rounded, rows_by, variant_of, write_rows, wrote

app = typer.Typer(add_completion=False)

# How many standard errors two means must be apart before the difference is called a difference.
# Two, so the claim is roughly 95% rather than "the means differ", which at these spreads they
# always do.
SIGMA = 2.0

# What a table's `drop_sem` is the standard error *of*, written into every row so that a
# consumer reading one CSV never has to know which command produced it.
OWN_BASIS = "own_samples"
COHORT_BASIS = "cohort_between_runs"

# The two conditions ranking tier 1 is the difference of. Written into one pass so the
# difference is paired within an invocation, exactly as the control's bar is.
TIER1_CONDITIONS = (WILDTYPE, COMPLEX)


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


def ranked_designs(dependence: Path, column: str = "drop_mean") -> list[tuple[str, float]]:
    """Designs by peptide-dependence drop, best first, with the control left out.

    The column is named because the campaign has two table shapes: the pilot screen wrote
    `peptide_dependence_drop` per design, and every pass since writes `drop_mean` over repeated
    samples. Defaulting to the latter means a caller that forgets gets the measurement of
    record, and the one retired caller says which it wants.
    """
    ranked = sorted(
        (
            (design, float(row[column]))
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
    ranked = ranked_designs(dependence, column="peptide_dependence_drop")
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


def difference_sem(high: Sequence[float], low: Sequence[float]) -> float:
    """Standard error of `mean(high) - mean(low)` for two independently sampled conditions.

    Each condition's mean carries `sd / sqrt(n)` and the difference is their difference, so the
    two combine in quadrature. The division happens per condition rather than once at the end
    because the conditions are sampled independently and may carry different `n`. Both the
    peptide-dependence drop and the tier-1 wild-type difference are this shape, so the rule
    lives here and neither restates it.
    """
    return math.hypot(
        stdev(high) / math.sqrt(len(high)),
        stdev(low) / math.sqrt(len(low)),
    )


@dataclass(frozen=True)
class Contrast:
    """One design in one invocation: the difference of two conditions' repeated draws.

    `Spread` is the same shape for the drop's own pair of conditions and carries that pair's
    vocabulary in its column names. This is the condition-agnostic form the tier-1 pass needs,
    so `Paired` and `Cohort` can supervise it without a second statistics path beside them.
    """

    design: str
    replicate: str
    high: tuple[float, ...]
    low: tuple[float, ...]

    @property
    def mean(self) -> float:
        return mean(self.high) - mean(self.low)

    @property
    def sem(self) -> float:
        return difference_sem(self.high, self.low)

    @property
    def samples(self) -> str:
        return f"{len(self.high)}/{len(self.low)}"


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

        This, not `drop_sd`, is the denominator of the tie rule. See `difference_sem`, which
        is the same rule the tier-1 pass uses.
        """
        return difference_sem(self.complexed, self.apo)

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
            "sem_basis": OWN_BASIS,
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
    dependence: Path = typer.Option(
        Path("cofold_spread_pooled.csv"),
        help="the measurement of record, whose column ranked_designs defaults to",
    ),
    targets: Path = typer.Option(Path("targets/cofold_targets.json")),
    gen_out: Path = typer.Option(Path("gen_out")),
    controls: Path = typer.Option(Path("targets/controls.json")),
    out: Path = typer.Option(Path("cofold_wt_specs")),
    count: int = typer.Option(20, help="finalists, which is the submission cap"),
) -> None:
    """Specs for ranking tier 1: both peptides, for the finalists and for the control.

    Both conditions in the same pass, and the control alongside them, because the difference
    tier 1 asks for is only interpretable paired. The control's own bar is `complex_wt` minus
    `complex_g12d` measured inside one invocation, so the invocation's offset cancels out of it;
    a design's difference has to be formed the same way or it is not comparable to the bar. An
    earlier version of this pass wrote `complex_wt` alone and took the G12D side from a different
    invocation's table. On the control the two conditions move together within a run (+0.0024,
    +0.0023) while the run itself shifts both by ~0.0012 -- 44% of the bar -- so that difference
    carried a run offset the bar does not contain, and the error bar did not know about it.
    """
    chosen = ranked_designs(dependence)[:count]
    finalists = [name for name, _ in chosen]
    conditions = conditions_from(targets)
    for needed in TIER1_CONDITIONS:
        if needed not in conditions:
            raise typer.BadParameter(f"{targets} has no condition {needed}")

    tcr = json.loads(controls.read_text())["tcr"]
    written = write_entries(
        out,
        [
            *design_entries(finalists, design_sequences(gen_out), conditions, TIER1_CONDITIONS),
            *control_entries(tcr, conditions, TIER1_CONDITIONS),
        ],
    )
    typer.echo(
        f"{len(chosen)} finalists by drop,"
        f" {chosen[0][1]:+.3f} down to {chosen[-1][1]:+.3f}"
    )
    typer.echo(
        f"wrote {written} specs to {out}:"
        f" {len(finalists)} designs plus {CONTROL} through {', '.join(TIER1_CONDITIONS)}."
        " Run it at least twice, as separate invocations -- a single run cannot see the"
        " between-run spread of a difference of two sampled conditions"
    )


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


@dataclass(frozen=True)
class Cohort:
    """The whole field's reproducibility, which is the denominator each design has to borrow.

    Two invocations cannot estimate a design's own between-run spread: the estimate carries
    about 70% relative error, and the pair that happens to agree reports zero. Taken at face
    value that makes a +0.011 drop look significant and a reproducible +0.112 look infinitely
    so. Pooling across designs buys roughly one degree of freedom each, and every per-design
    verdict below divides by the pooled number rather than by its own two points.
    """

    between_sd: float
    within_sem: float

    @classmethod
    def of(
        cls, within_sems: Mapping[str, Sequence[float]], paired: Mapping[str, Paired]
    ) -> Cohort:
        """`within_sems` is each design's per-invocation within-run error, whatever produced it.

        Taking the numbers rather than the rows is what lets the drop and the tier-1 difference
        share one denominator rule instead of each growing its own.
        """
        designs = [design for design in paired if design != CONTROL] or list(paired)
        return cls(
            between_sd=pooled_between_sd([paired[design].between_sd for design in designs]),
            within_sem=mean(sem for design in designs for sem in within_sems[design]),
        )

    @property
    def run_effect(self) -> float:
        return run_effect_sd(self.between_sd, self.within_sem)

    def sem(self, runs: int) -> float:
        """What a mean of `runs` invocations is worth, on the cohort's measured spread."""
        return self.between_sd / math.sqrt(runs)

    def verdict(self, drops: Paired, against: float = 0.0, against_sem: float = 0.0) -> str:
        return comparison(drops.mean, self.sem(drops.runs), against, against_sem)


def echo_replicate_table(
    by_design: Mapping[str, list[RunRow]], paired: Mapping[str, Paired], cohort: Cohort
) -> None:
    typer.echo(
        f"\n  {'design':20s} {'runs':>4s} {'drop':>8s} {'own sd':>7s}"
        f" {'within':>7s} {'sem':>6s} {'vs 0':>6s}  per-run drops"
    )
    for design, drops in sorted(paired.items(), key=lambda kv: -kv[1].mean):
        within = mean(row.spread.drop_sem for row in by_design[design])
        typer.echo(
            f"  {design:20s} {drops.runs:4d} {drops.mean:+8.3f} {drops.between_sd:7.3f}"
            f" {within:7.3f} {cohort.sem(drops.runs):6.3f}"
            f" {cohort.verdict(drops):>6s}"
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
    found = Replicated.of(out_dir, "a between-run spread needs two independent invocations")
    walk, by_design = found.walk, found.by_design
    typer.echo(
        f"{walk.directories} prediction directories over {found.invocations} invocations:"
        f" {len(walk.rows)} design-runs scored, {walk.unscored} with no samples,"
        f" {walk.unnamed} unrecognised names"
    )
    singles = len({row.spread.design for row in walk.rows}) - len(by_design)
    typer.echo(
        f"{len(by_design)} designs measured more than once;"
        f" {singles} appear in one invocation only and are excluded"
    )

    paired, cohort = found.paired, found.cohort
    echo_replicate_table(by_design, paired, cohort)
    echo_bar(by_design, paired, cohort)
    echo_wildtype(walk.rows)
    typer.echo(wrote(out, write_rows(walk.rows, out)))


def echo_bar(
    by_design: Mapping[str, list[RunRow]], paired: Mapping[str, Paired], cohort: Cohort
) -> None:
    """The control's drop, and how the field stands against it on the cohort's own spread."""
    typer.echo(
        f"\ncohort reproducibility over {len(paired)} designs: between-run SD"
        f" {cohort.between_sd:.3f} against a within-run error of {cohort.within_sem:.3f},"
        f" leaving a run effect of {cohort.run_effect:.3f}"
    )
    if (control := paired.get(CONTROL)) is None:
        typer.echo(
            f"NO CONTROL: {CONTROL} is not repeated here, so this pass calibrates nothing."
            " Nothing downstream may rank against it; re-run with the control's specs included."
        )
        return
    within = mean(row.spread.drop_sem for row in by_design[CONTROL])
    typer.echo(
        f"BAR: the control's drop is {control.mean:+.3f} over {control.runs} invocations,"
        f" spread {control.between_sd:.3f} against a within-run error of {within:.3f}"
    )

    designs = [p for design, p in paired.items() if design != CONTROL]
    if not designs:
        return
    bar_sem = cohort.sem(control.runs)
    verdicts = [cohort.verdict(p, control.mean, bar_sem) for p in designs]
    tally = {v: verdicts.count(v) for v in ("above", "tied", "below")}
    typer.echo(
        f"against it, on the pooled spread: {tally['above']} above, {tally['tied']} tied,"
        f" {tally['below']} below -- and {sum(cohort.verdict(p) == 'above' for p in designs)}"
        f" of {len(designs)} clear zero"
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


@dataclass(frozen=True)
class Replicated:
    """Every design measured by more than one invocation, with the denominator they share.

    Both replicate commands open with the same five steps and the same refusal, so they are
    established once here: a change to how the cohort is pooled cannot then apply to one
    command and not the other.
    """

    walk: Walk
    by_design: dict[str, list[RunRow]]
    paired: dict[str, Paired]
    cohort: Cohort

    @classmethod
    def of(cls, out_dirs: Sequence[Path], advice: str) -> Replicated:
        walk = run_rows(out_dirs)
        by_design = repeated(walk.rows)
        if not by_design:
            raise typer.BadParameter(
                f"no design appears in more than one of {[str(d) for d in out_dirs]}; {advice}"
            )
        paired = paired_drops(by_design)
        return cls(
            walk=walk,
            by_design=by_design,
            paired=paired,
            cohort=Cohort.of(
                {d: [row.spread.drop_sem for row in runs] for d, runs in by_design.items()},
                paired,
            ),
        )

    @property
    def invocations(self) -> int:
        return len({row.replicate for row in self.walk.rows})


@dataclass(frozen=True)
class Pooled:
    """One design's drop over every invocation, in the shape `report-spread` writes.

    `drop_sd` keeps that table's column name because `submission.candidates` reads it, and
    keeping it means `drop_sem == drop_sd / sqrt(n)` still holds -- but it holds the cohort's
    between-run spread, so it is **the same constant on every row** here. `sem_basis` says so in
    the data rather than leaving a reader to infer it from which command wrote the file.

    `drop_own_sd` is the design's own spread across invocations, carried for audit and never
    divided by: at two runs it is what `Cohort` exists to stop anyone trusting. It is also the
    headline for a design that reproduces -- 0.0001 on the leader -- which is why it travels as
    far as the evidence file instead of stopping at `cofold_replicates.csv`.

    The control row gets the designs' pooled scale like every other row, not its own 10-run
    spread of 0.016. That spread is inflated by the one run whose apo pose degraded, and using
    it would make the table's error bars incommensurable between rows; it moves the margin by
    1.4% and flips no verdict.
    """

    design: str
    drops: Paired
    cohort: Cohort

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "runs": self.drops.runs,
            "drop_mean": round(self.drops.mean, 4),
            "drop_sd": round(self.cohort.between_sd, 4),
            "drop_sem": round(self.cohort.sem(self.drops.runs), 4),
            "sem_basis": COHORT_BASIS,
            "drop_significant": self.cohort.verdict(self.drops) == "above",
            "drop_own_sd": round(self.drops.between_sd, 4),
            "drop_min": round(min(self.drops.values), 4),
            "drop_max": round(max(self.drops.values), 4),
        }


@app.command()
def pool_replicates(
    out_dir: list[Path] = typer.Option(
        [Path("cofold_full5_out")],
        help="repeatable: every pass whose runs should be pooled into one ranking",
    ),
    out: Path = typer.Option(Path("cofold_spread_pooled.csv")),
) -> None:
    """Pool repeated invocations into the table the submission ranks on."""
    found = Replicated.of(out_dir, "rank on report-spread's single-pass table instead")
    rows = sorted(
        (Pooled(design, drops, found.cohort) for design, drops in found.paired.items()),
        key=lambda pooled: -pooled.drops.mean,
    )
    typer.echo(
        f"pooled {len(rows)} designs over {found.invocations} invocations;"
        f" every drop_sem is the cohort's {found.cohort.between_sd:.3f} spread over sqrt(runs)"
    )
    typer.echo(wrote(out, write_rows(rows, out)))


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


def tier1_contrasts(out_dirs: Sequence[Path]) -> tuple[list[Contrast], Walk]:
    """Every design's wild-type-minus-G12D difference, once per independent invocation.

    Paired inside the invocation that produced both conditions, which is the whole point of the
    pass: an invocation shifts both conditions together, so forming the difference within it
    cancels that shift. Taking the two sides from different invocations does not, and the error
    bar cannot see what it left out.
    """
    contrasts: list[Contrast] = []
    directories = unscored = unnamed = 0
    for out_dir in out_dirs:
        for predictions in invocation_predictions(out_dir):
            screen = collect(predictions, sample_iptms)
            directories += screen.directories
            unscored += screen.unscored
            unnamed += screen.unnamed
            for design, conditions in sorted(screen.runs.items()):
                high, low = conditions.get(WILDTYPE, ()), conditions.get(COMPLEX, ())
                if len(high) > 1 and len(low) > 1:
                    contrasts.append(
                        Contrast(
                            design=design,
                            replicate=predictions.parent.parent.name,
                            high=high,
                            low=low,
                        )
                    )
    walk = Walk(rows=(), directories=directories, unscored=unscored, unnamed=unnamed)
    return contrasts, walk


@dataclass(frozen=True)
class WildtypeRow:
    """One design's tier-1 difference against the control's, pooled over invocations.

    The bar is the control measured in the same invocations, so nothing is carried in from
    another pass -- see the module docstring for why only the comparison against it is readable.
    """

    design: str
    differences: Paired
    cohort: Cohort
    bar: Paired

    @property
    def sem(self) -> float:
        """The cohort's pooled between-run spread over this design's runs, as the drop uses."""
        return self.cohort.sem(self.differences.runs)

    @property
    def bar_sem(self) -> float:
        return self.cohort.sem(self.bar.runs)

    @property
    def versus_bar(self) -> str:
        return self.cohort.verdict(self.differences, self.bar.mean, self.bar_sem)

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "variant": variant_of(self.design),
            "runs": self.differences.runs,
            "wt_minus_g12d": rounded(self.differences.mean),
            "wt_minus_g12d_sd": rounded(self.cohort.between_sd),
            "wt_minus_g12d_sem": rounded(self.sem),
            "sem_basis": COHORT_BASIS,
            "wt_minus_g12d_own_sd": rounded(self.differences.between_sd),
            "control_bar": rounded(self.bar.mean),
            "control_bar_sem": rounded(self.bar_sem),
            "control_bar_runs": self.bar.runs,
            "versus_bar": self.versus_bar,
            "per_run": " ".join(f"{v:+.4f}" for v in self.differences.values),
        }


def echo_screen(label: str, directories: int, scored: int, unscored: int, unnamed: int) -> None:
    """The denominators every pass reports. One wording, so two commands cannot diverge."""
    typer.echo(
        f"{directories} prediction directories in {label}: {scored} scored,"
        f" {unscored} with no samples, {unnamed} unrecognised names"
    )


@app.command()
def report_wildtype(
    out_dir: list[Path] = typer.Option(
        [Path("cofold_wt_out")],
        help="repeatable: a directory of replicate sub-directories, or one pass's --out_dir",
    ),
    out: Path = typer.Option(Path("cofold_wildtype.csv")),
) -> None:
    """Ranking tier 1, paired within each invocation and scaled by the cohort. Never gated.

    The difference has to be formed inside each invocation and needs at least two invocations
    to scale it; this command requires both rather than annotating their absence.
    """
    contrasts, walk = tier1_contrasts(out_dir)
    echo_screen(
        ", ".join(str(d) for d in out_dir),
        walk.directories,
        len(contrasts),
        walk.unscored,
        walk.unnamed,
    )

    by_design: dict[str, list[Contrast]] = {}
    for contrast in contrasts:
        by_design.setdefault(contrast.design, []).append(contrast)
    by_design = {d: runs for d, runs in by_design.items() if len(runs) > 1}
    if CONTROL not in by_design:
        raise typer.BadParameter(
            f"{CONTROL} has no repeated {WILDTYPE}/{COMPLEX} pair in"
            f" {[str(d) for d in out_dir]}, so these differences have no bar. Tier 1 is only"
            " interpretable against the one receptor known to discriminate these peptides;"
            " re-run the pass with select-wildtype, which writes the control's specs beside the"
            " finalists'."
        )

    paired = {
        design: Paired(design, tuple(c.mean for c in runs)) for design, runs in by_design.items()
    }
    cohort = Cohort.of({d: [c.sem for c in runs] for d, runs in by_design.items()}, paired)
    bar = paired[CONTROL]
    rows = [
        WildtypeRow(design=design, differences=differences, cohort=cohort, bar=bar)
        for design, differences in paired.items()
        if design != CONTROL
    ]

    typer.echo(
        f"\nBAR: the control's own {WILDTYPE} minus {COMPLEX} is {bar.mean:+.4f}"
        f" +- {cohort.sem(bar.runs):.4f} over {bar.runs} invocations"
        f" ({' '.join(f'{v:+.4f}' for v in bar.values)}) -- {bar.verdict} zero."
        " Positive means this predictor prefers the peptide the control does NOT recognise,"
        " so a design's position relative to this bar is the only readable thing here"
    )
    typer.echo(
        f"SCALE: cohort between-run SD {cohort.between_sd:.4f} against a within-run"
        f" {cohort.within_sem:.4f}, leaving a run effect of {cohort.run_effect:.4f}"
    )

    typer.echo(
        f"\n  {'design':20s} {'variant':14s} {'runs':>4s} {'wt-G12D':>9s} {'own sd':>7s}"
        f" {'sem':>7s} {'vs bar':>7s}  per-run"
    )
    for row in sorted(rows, key=lambda r: r.differences.mean):
        typer.echo(
            f"  {row.design:20s} {variant_of(row.design):14s} {row.differences.runs:4d}"
            f" {row.differences.mean:+9.4f} {row.differences.between_sd:7.4f} {row.sem:7.4f}"
            f" {row.versus_bar:>7s}  {' '.join(f'{v:+.4f}' for v in row.differences.values)}"
        )

    verdicts = [row.versus_bar for row in rows]
    tally = {v: verdicts.count(v) for v in ("above", "tied", "below")}
    across = [row.differences.mean for row in rows]
    typer.echo(
        f"\n{tally['below']} of {len(rows)} finalists sit below the bar, {tally['tied']} tied,"
        f" {tally['above']} above. Across designs -- a spread over designs, not over runs --"
        f" the difference averages {mean(across):+.4f}"
        f" +- {stdev(across) / math.sqrt(len(across)):.4f}"
    )
    typer.echo(
        "REPORTED, NOT GATED: the bar's sign is inverted on the one molecule that can calibrate"
        " it, and one predictor's agreement with itself across runs is not a second method."
        " No design is ranked on this column."
    )
    typer.echo(wrote(out, write_rows(rows, out)))


if __name__ == "__main__":
    app()
