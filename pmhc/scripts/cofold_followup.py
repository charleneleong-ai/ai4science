"""The co-fold passes that run after the screen has ranked the survivors.

**Sampling noise.** The screen draws one diffusion sample per condition, so a design's drop is
a difference of two single draws and carries the variance of both. Ranking 78 designs on it is
only defensible once the spread of that difference is known, and the screen cannot report it:
one sample has no spread. `select-noise` writes specs for a subset spanning the observed drop
ranking, to be run at several samples, and `report-noise` gives each design the standard
deviation of its per-condition ipTM plus the implied spread of the drop. Designs whose drops
differ by less than that spread are ties, and are reported as ties.

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
def report_noise(
    out_dir: Path = typer.Option(Path("cofold_noise_out"), help="boltz --out_dir"),
    out: Path = typer.Option(Path("cofold_noise.csv")),
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

    typer.echo(
        f"\n  {'design':20s} {'n':>5s} {'complex':>15s} {'apo':>15s}"
        f" {'drop':>8s} {'+-':>6s} {'range':>16s}"
    )
    for row in sorted(rows, key=lambda r: -r.drop_mean):
        typer.echo(
            f"  {row.design:20s} {row.samples:>5s}"
            f" {row.complex_mean:7.3f}+-{row.complex_sd:<6.3f}"
            f" {row.apo_mean:7.3f}+-{row.apo_sd:<6.3f}"
            f" {row.drop_mean:+8.3f} {row.drop_sd:6.3f}"
            f" [{row.drop_min:+.3f},{row.drop_max:+.3f}]"
        )

    worst = max(r.drop_sd for r in rows)
    control = next((r for r in rows if r.design == CONTROL), None)
    typer.echo(
        f"\nlargest drop spread over {len(rows)} entities: {worst:.3f}."
        " Two designs whose drops differ by less than this are tied."
    )
    if control is not None:
        typer.echo(
            f"the control's own drop is {control.drop_mean:+.3f} +- {control.drop_sd:.3f}"
            f" over {control.samples} samples, so its bar is a band, not a line"
        )
    typer.echo(wrote(out, write_rows(rows, out)))


if __name__ == "__main__":
    app()
