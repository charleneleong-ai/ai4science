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
reported, never gated: the control measures 0.957 on G12D against 0.958 on wild type, so the one
receptor known to discriminate these peptides does not discriminate them here, and no difference
this predictor reports on a design could mean anything. The column exists so that the claim is
auditable rather than asserted.
"""

from __future__ import annotations

import json
import math
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
    predictions_dir,
    write_entries,
)
from common import rows_by, write_rows, wrote

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
