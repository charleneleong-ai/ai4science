"""Render the three charts the PR and the writeup are read from.

One chart per question, each drawn from the measured tables rather than from numbers retyped
into the plotting code:

1. **Why burial and not carboxylates.** The linkage against distance to the target, both pKa
   methods, over every position measured on both the de novo pool and the 3DS18 control. This is
   the campaign's central claim and the only chart that argues a mechanism.
2. **Where the twenty slots went.** Each chosen design's switch, grouped by chassis, against the
   two log units that "no detectable binding at pH 6.0" would need.
3. **The funnel, with every denominator.** 300 generated to 20 submitted, so the three
   survived-of-N statements in the writeup can be checked against one picture.

Rendered on an explicit light surface rather than a transparent one: a committed PNG cannot
adapt, and transparent backgrounds put dark axis text on GitHub's dark theme. Colours are the
first two categorical slots of the validated default palette, which clear the all-pairs CVD and
normal-vision floors needed for a scatter (protan dE 24.7, normal dE 33.6).
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated

import matplotlib
import typer

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be set before pyplot loads)
from csvio import number, read_rows  # noqa: E402

SURFACE = "#fcfcfb"
INK, INK_SOFT, GRID = "#0b0b0b", "#52514e", "#dcdbd6"
PROPKA, PKAI = "#2a78d6", "#eb6834"
SINGLE = "#2a78d6"

# What "no detectable binding at pH 6.0" costs, from `ph_switch.release_log_units`.
TARGET_LOG_UNITS = 2.0

app = typer.Typer(add_completion=False)


def style(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SOFT, labelsize=9)


def figure(width: float, height: float) -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=(width, height), facecolor=SURFACE)
    style(ax)
    return fig, ax


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    typer.echo(f"wrote {path}")


def short_chassis(parent: str) -> str:
    """The chassis name without the spec prefix the pooled table doubles."""
    return (parent.replace("conserved_s_conserved_s_", "cons_s_")
            .replace("conserved_conserved_", "cons_")
            .replace("carboxylate_s_carboxylate_s_", "carb_s_"))


def carboxylate_distance(cell: str) -> float | None:
    """The angstroms out of a `GLU74@5.6A` cell, or None when no acid was within range."""
    if not cell or "@" not in cell:
        return None
    return number(cell.rsplit("@", 1)[1].removesuffix("A"), 0.0) or None


def plot_burial(rows: Sequence[dict[str, str]], floor: float, out: Path) -> None:
    """The same linkage against each candidate explanation, so the contrast is the chart.

    Two panels rather than one, because the claim is comparative: a single panel against
    distance-to-target shows a relationship but cannot show that the carboxylate story fails,
    and a title asserting it while plotting only burial would be arguing past the picture.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), facecolor=SURFACE, sharey=True)
    panels = (
        (axes[0], "distance to the target (Å)", lambda r: number(r["target_distance"])
         if r.get("target_distance") else None),
        (axes[1], "distance to the nearest target carboxylate (Å)",
         lambda r: carboxylate_distance(r.get("nearest_carboxylate", ""))),
    )

    for ax, xlabel, distance in panels:
        style(ax)
        for method, colour, label in (("propka", PROPKA, "PROPKA 3.5.1"),
                                      ("pkai", PKAI, "pKAI 1.2.0")):
            points = [(d, number(r[f"delta_{method}"])) for r in rows
                      if (d := distance(r)) is not None and r.get(f"delta_{method}")]
            ax.scatter([x for x, _ in points], [y for _, y in points], s=22, c=colour,
                       edgecolors=SURFACE, linewidths=0.7, alpha=0.8, label=label, zorder=3)
        ax.axhline(0, color=INK_SOFT, linewidth=1.0, zorder=2)
        ax.axhline(-floor, color=INK_SOFT, linewidth=1.2, linestyle="--", zorder=2)
        ax.set_xlabel(xlabel, color=INK, fontsize=10)

    axes[0].set_ylabel("Δ = pKa(complex) − pKa(free)", color=INK, fontsize=10)
    axes[0].annotate(f"coupling floor −{floor:g}", xy=(axes[0].get_xlim()[0], -floor),
                     xytext=(6, -13), textcoords="offset points", ha="left", fontsize=8.5,
                     color=INK_SOFT)
    axes[1].annotate("above 0: binds harder at pH 6.0", xy=(axes[1].get_xlim()[1], 0),
                     xytext=(-6, 6), textcoords="offset points", ha="right", fontsize=8.5,
                     color=INK_SOFT)
    axes[0].annotate("native His69,\nthe internal control", xy=(11.0, -0.02),
                     xytext=(-10, -46), textcoords="offset points", ha="right", fontsize=8.5,
                     color=INK_SOFT,
                     arrowprops={"arrowstyle": "-", "color": INK_SOFT, "linewidth": 0.8})
    axes[1].legend(frameon=False, fontsize=9, labelcolor=INK_SOFT, loc="lower right")

    fig.suptitle("Burial sets the switch; proximity to a carboxylate does not",
                 color=INK, fontsize=12.5, x=0.065, ha="left", y=1.0)
    save(fig, out)


def plot_chosen(rows: Sequence[dict[str, str]], out: Path) -> None:
    """Each chosen design's switch, blocked by chassis so the spread of slots is visible."""
    ordered = sorted(rows, key=lambda r: (r["parent"], -number(r["switch_log_units"])))

    # One blank row between chassis, so each group's name has somewhere to sit that is not on
    # top of a bar.
    positions, labels, values, headings = [], [], [], []
    cursor, previous = 0.0, None
    for row in ordered:
        if previous is not None and row["parent"] != previous:
            cursor += 1.0
        if row["parent"] != previous:
            headings.append((short_chassis(row["parent"]), cursor - 0.85))
            previous = row["parent"]
        positions.append(cursor)
        labels.append(row["name"].split("-", 1)[1])
        values.append(number(row["switch_log_units"]))
        cursor += 1.0

    fig, ax = figure(8.4, 7.0)
    ax.barh(positions, values, height=0.68, color=SINGLE, zorder=3)
    for y, value in zip(positions, values, strict=True):
        fold = 10 ** value
        ax.annotate(f"{fold:.0f}×" if fold >= 10 else f"{fold:.1f}×",
                    xy=(value, y), xytext=(4, 0), textcoords="offset points", va="center",
                    fontsize=8.5, color=INK_SOFT)
    for name, y in headings:
        ax.annotate(name, xy=(0, y), xytext=(1, 0), textcoords="offset points", fontsize=9,
                    color=INK, va="center", fontweight="bold")

    ax.axvline(TARGET_LOG_UNITS, color=INK_SOFT, linewidth=1.2, linestyle="--", zorder=2)
    ax.annotate("2 log units: no detectable binding at pH 6.0",
                xy=(TARGET_LOG_UNITS, cursor / 2), xytext=(-5, 0), textcoords="offset points",
                ha="center", va="center", fontsize=8.5, color=INK_SOFT, rotation=90)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(cursor - 0.4, -1.4)
    ax.set_xlim(0, TARGET_LOG_UNITS * 1.16)
    ax.set_xlabel("predicted switch (log₁₀ fold weaker at pH 6.0)", color=INK,
                  fontsize=10)
    ax.set_title(f"The {len(ordered)} submitted designs, by chassis",
                 color=INK, fontsize=12, loc="left", pad=10)
    save(fig, out)


Stages = Sequence[tuple[str, int, int]]


def plot_funnel(generated: Stages, engineered: Stages, out: Path) -> None:
    """Every stage with its denominator, so no survived-of-N claim rests on prose alone.

    Two panels because there are two populations and one bar chart cannot hold both: the
    generator's output is counted in designs, and what is submitted is counted in histidine
    variants *of* the designs that survived. Plotting them on one axis makes the third bar rise
    above the second, which reads as a funnel that gains candidates.
    """
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 5.0), facecolor=SURFACE,
                             gridspec_kw={"height_ratios": [3, 3]})
    for ax, stages, unit in ((axes[0], generated, "designs the generator produced"),
                             (axes[1], engineered, "histidine variants of the survivors")):
        style(ax)
        positions = list(range(len(stages)))
        ax.barh(positions, [kept for _, kept, _ in stages], height=0.6, color=SINGLE, zorder=3)
        for y, (_, kept, total) in zip(positions, stages, strict=True):
            ax.annotate(f"{kept} of {total}", xy=(kept, y), xytext=(5, 0),
                        textcoords="offset points", va="center", fontsize=9, color=INK_SOFT)
        ax.set_yticks(positions)
        ax.set_yticklabels([label for label, _, _ in stages], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlim(0, max(total for _, _, total in stages) * 1.26)
        ax.set_xlabel(unit, color=INK, fontsize=9.5)

    fig.suptitle("The funnel, with every denominator", color=INK, fontsize=12.5, x=0.125,
                 ha="left", y=1.0)
    fig.tight_layout(h_pad=2.2)
    save(fig, out)


@app.command()
def main(
    switch: Annotated[Path, typer.Option(help="de novo per design/position agreement CSV")],
    control: Annotated[Path, typer.Option(help="3DS18 control agreement CSV")],
    chosen: Annotated[Path, typer.Option(help="the selected designs")],
    rejected: Annotated[Path, typer.Option(help="every candidate with its reason")],
    triage: Annotated[Path, typer.Option(help="per-design triage verdicts")],
    thresholds: Annotated[Path, typer.Option(help="the verifier's resolved thresholds")],
    out_dir: Annotated[Path, typer.Option(help="where the PNGs are written")],
) -> None:
    floor = float(json.loads(thresholds.read_text())["coupling_floor"])
    linkage = read_rows(switch) + read_rows(control)
    candidates = read_rows(rejected)
    triaged = read_rows(triage)
    picks = read_rows(chosen)

    chassis = {r["parent"] for r in candidates}
    admitted = sum(not r["reason"] or r["reason"] == "not in the top slots"
                   for r in candidates)

    plot_burial(linkage, floor, out_dir / "burial-vs-linkage.png")
    plot_chosen(picks, out_dir / "chosen-designs.png")
    plot_funnel(
        [("generated by BoltzGen", len(triaged), len(triaged)),
         ("pass composition triage", sum(r["passes"] == "yes" for r in triaged), len(triaged)),
         ("chosen as a chassis", len(chassis), len(triaged))],
        [("built and measured", len(candidates), len(candidates)),
         ("admissible after the gates", admitted, len(candidates)),
         ("submitted", sum(r["chosen"] == "True" for r in candidates), len(candidates))],
        out_dir / "funnel.png")


if __name__ == "__main__":
    app()
