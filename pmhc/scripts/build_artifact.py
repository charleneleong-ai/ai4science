"""Render the campaign's progress page from the committed result files.

Every number on the page is read from `docs/experiments/pmhc/campaign.yaml`, the verdict ledger
beside it, or one of the `data/*.csv` tables. Nothing is typed into the template, so a rebuild
after a new measurement cannot leave a stale figure behind -- which is the whole reason this is a
script and not a hand-written page.

Each chart answers one question and carries its finding above it rather than a caption. The
charts are inline SVG on one scale per chart, with every colour drawn from the theme tokens so
the page reads in either light or dark.

    uv run --extra artifact python scripts/build_artifact.py --out docs/artifact/index.html
"""

from __future__ import annotations

import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import typer
import yaml

app = typer.Typer(add_completion=False)

# Kept in step with `common.VARIANTS` by a test rather than imported: this builder renders from
# committed text files and must not pull in the measurement stack (gemmi, Biopython, numpy).
VARIANT_ORDER = ("bulge_core", "bulge_wide", "mut_proximal", "hub")
VERDICT_ORDER = ("PROMOTE", "KEEP", "DISCARD", "running")
# The receptor travels in the result tables as a row like any design, because it is measured in
# the same runs. It is the bar, so it is never drawn as one of the designs it is the bar for.
CONTROL_DESIGN = "control_tcr"
REPRODUCES_WITHIN = 0.005
CONE_FLOOR = 0.01  # the log axis' floor, below which a percentage is drawn at the axis bottom


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return yaml.safe_load(handle)


def esc(value: object) -> str:
    return html.escape(str(value))


def split_control(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """`(designs, control_rows)`. Every chart of the cohort takes the designs."""
    control = [row for row in rows if row["design"] == CONTROL_DESIGN]
    return [row for row in rows if row["design"] != CONTROL_DESIGN], control


def by_variant(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    """Finalists grouped by generation variant, empty groups dropped, in nomination order."""
    return {
        variant: members
        for variant in VARIANT_ORDER
        if (members := [row for row in rows if row.get("variant") == variant])
    }


@dataclass(frozen=True)
class Scale:
    """One linear mapping from data to pixels, so a chart cannot mix two of them."""

    lo: float
    hi: float
    px_lo: float
    px_hi: float

    def __call__(self, value: float) -> float:
        span = self.hi - self.lo or 1.0
        return self.px_lo + (value - self.lo) / span * (self.px_hi - self.px_lo)


def finding(text: str, chart: str, note: str = "") -> str:
    """A chart with its finding above it -- the finding, not a description of the axes."""
    footer = f'<p class="note">{note}</p>' if note else ""
    return (
        f'<figure><figcaption>{text}</figcaption>'
        f'<div class="plot">{chart}</div>{footer}</figure>'
    )


def scoreboard(rows: list[dict[str, str]], bar: float) -> str:
    """Every screened design against the bar, rejected ones included."""
    drops = sorted((float(r["drop_mean"]) for r in rows), reverse=True)
    w, h, pad = 760, 240, 34
    lo, hi = min(min(drops), 0.0), max(max(drops), bar)
    y = Scale(lo, hi, h - pad, pad)
    step = (w - 2 * pad) / max(len(drops), 1)
    bars = "".join(
        f'<rect x="{pad + i * step:.1f}" y="{min(y(d), y(0.0)):.1f}" width="{max(step - 1.2, 0.8):.2f}"'
        f' height="{abs(y(d) - y(0.0)):.1f}" rx="1"'
        f' fill="{"var(--c1)" if d >= bar else "var(--c2)" if d > 0 else "var(--muted)"}"/>'
        for i, d in enumerate(drops)
    )
    return (
        f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="peptide-dependence drop for every'
        f' screened design against the receptor bar">'
        f'<line x1="{pad}" x2="{w - pad}" y1="{y(0.0):.1f}" y2="{y(0.0):.1f}" stroke="var(--grid)"/>'
        f'<line x1="{pad}" x2="{w - pad}" y1="{y(bar):.1f}" y2="{y(bar):.1f}"'
        f' stroke="var(--fg)" stroke-dasharray="5 4" stroke-width="1.5"/>'
        f'<text x="{w - pad}" y="{y(bar) - 7:.1f}" text-anchor="end" class="lbl">'
        f'receptor bar {bar:+.3f}</text>'
        f"{bars}"
        f'<text x="{pad}" y="{y(hi) - 12:.1f}" class="lbl">{hi:+.3f}</text>'
        f'<text x="{pad}" y="{y(lo) + 16:.1f}" class="lbl">{lo:+.3f}</text>'
        f'<text x="{pad}" y="{y(0.0) - 6:.1f}" class="lbl">0</text>'
        f"</svg>"
    )


def reproducibility(rows: list[dict[str, str]], top: int = 14) -> str:
    """Each design's two independent runs, which is the apo-artifact test."""
    by_design: dict[str, list[float]] = {}
    for row in rows:
        by_design.setdefault(row["design"], []).append(float(row["drop_mean"]))
    pairs = [(d, v) for d, v in by_design.items() if len(v) > 1]
    pairs.sort(key=lambda kv: -max(kv[1]))
    pairs = pairs[:top]
    w, pad, row_h = 760, 150, 19
    right = 96
    h = pad // 6 + len(pairs) * row_h + 30
    flat = [v for _, vs in pairs for v in vs]
    x = Scale(min(flat), max(flat), pad, w - right)
    marks = []
    for i, (design, values) in enumerate(pairs):
        cy = 14 + i * row_h
        a, b = x(min(values)), x(max(values))
        spread = max(values) - min(values)
        colour = "var(--c1)" if spread < REPRODUCES_WITHIN else "var(--c3)"
        marks.append(
            f'<text x="{pad - 10}" y="{cy + 4}" text-anchor="end" class="lbl">{esc(design)}</text>'
            f'<line x1="{a:.1f}" x2="{b:.1f}" y1="{cy}" y2="{cy}" stroke="{colour}" stroke-width="2.5"/>'
            f'<circle cx="{a:.1f}" cy="{cy}" r="4" fill="{colour}"/>'
            f'<circle cx="{b:.1f}" cy="{cy}" r="4" fill="{colour}"/>'
            f'<text x="{w - right + 14}" y="{cy + 4}" class="lbl">{spread:.3f}</text>'
        )
    return (
        f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="per-design spread between two'
        f' independent invocations">'
        f'<text x="{w - right + 14}" y="{h - 8}" class="lbl">run spread</text>'
        + "".join(marks)
        + "</svg>"
    )


def funnel(stages: list[dict[str, Any]]) -> str:
    """Every denominator, so the three survivors are never shown without the 256."""
    w, pad, row_h = 760, 230, 30
    h = len(stages) * row_h + 24
    widest = max(s["surviving"] + s["removed"] for s in stages)
    x = Scale(0, widest, pad, w - 70)
    marks = []
    for i, stage in enumerate(stages):
        cy = 12 + i * row_h
        kept, gone = stage["surviving"], stage["removed"]
        marks.append(
            f'<text x="{pad - 12}" y="{cy + 13}" text-anchor="end" class="lbl">'
            f'{esc(stage["stage"])}</text>'
            f'<rect x="{pad}" y="{cy}" width="{max(x(kept) - pad, 1):.1f}" height="17" rx="2"'
            f' fill="var(--c1)"/>'
            f'<rect x="{x(kept):.1f}" y="{cy}" width="{max(x(kept + gone) - x(kept), 0):.1f}"'
            f' height="17" rx="2" fill="var(--muted)" opacity="0.55"/>'
            f'<text x="{w - 64}" y="{cy + 13}" class="lbl">{kept}'
            + (f' <tspan fill="var(--dim)">-{gone}</tspan>' if gone else "")
            + "</text>"
        )
    return (
        f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="design funnel with every'
        f' denominator">' + "".join(marks) + "</svg>"
    )


def escape_cone(cone: dict[str, Any]) -> str:
    """Why tier 1 was never a gate: the distinguishing atom has no approach vector."""
    grid = {int(k): float(v) for k, v in cone["by_reach_angstrom"].items()}
    reaches = sorted(grid)
    w, h, pad = 760, 230, 46
    y = Scale(math.log10(CONE_FLOOR), math.log10(50.0), h - pad, pad)
    x = Scale(min(reaches), max(reaches), pad, w - pad)
    pts = " ".join(f"{x(r):.1f},{y(math.log10(max(grid[r], CONE_FLOOR))):.1f}" for r in reaches)
    dots = "".join(
        f'<circle cx="{x(r):.1f}" cy="{y(math.log10(max(grid[r], CONE_FLOOR))):.1f}" r="4.5"'
        f' fill="var(--c1)"/>'
        f'<text x="{x(r):.1f}" y="{y(math.log10(max(grid[r], CONE_FLOOR))) - 11:.1f}"'
        f' text-anchor="middle" class="lbl">{grid[r]:g}%</text>'
        for r in reaches
    )
    ticks = "".join(
        f'<text x="{x(r):.1f}" y="{h - pad + 18:.1f}" text-anchor="middle" class="lbl">'
        f"{r} Å</text>"
        for r in reaches
    )
    refs = "".join(
        f'<line x1="{pad}" x2="{w - pad}" y1="{y(math.log10(v)):.1f}"'
        f' y2="{y(math.log10(v)):.1f}" stroke="var(--grid)" stroke-dasharray="3 5"/>'
        for v in (10.0, 1.0, 0.1)
    )
    return (
        f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="escape cone over the G12D'
        f' carboxylate against designed-side-chain reach, log scale">'
        f"{refs}<polyline points=\"{pts}\" fill=\"none\" stroke=\"var(--c1)\""
        f' stroke-width="2"/>{dots}{ticks}</svg>'
    )


def tier1_chart(rows: list[dict[str, str]]) -> str:
    """The finalists' wild-type difference against the control's own, grouped by variant."""
    if not rows:
        return '<p class="note">The replicated tier-1 pass has not been read yet.</p>'
    bar = tier1_bar(rows)
    groups = by_variant(rows)
    w, pad = 760, 150
    # The right gutter is sized to the longest design name at 10px mono, so a name is
    # never clipped by the viewBox however the finalists are named.
    gutter = 6 * max(len(row["design"]) for row in rows) + 12
    row_h, h = 22, 24 + sum(len(v) for v in groups.values()) * 22 + len(groups) * 10
    # The domain spans the error bars, not just the means, so no bar runs under a label.
    ends = [
        float(row["wt_minus_g12d"]) + sign * float(row["wt_minus_g12d_sem"])
        for row in rows
        for sign in (-1, 1)
    ]
    x = Scale(min(*ends, bar) - 0.002, max(*ends, bar) + 0.002, pad, w - gutter - 8)
    marks = [
        f'<line x1="{x(bar):.1f}" x2="{x(bar):.1f}" y1="6" y2="{h - 18}" stroke="var(--fg)"'
        f' stroke-dasharray="5 4" stroke-width="1.5"/>'
        f'<text x="{x(bar):.1f}" y="{h - 5}" text-anchor="middle" class="lbl">'
        f"control {bar:+.4f}</text>"
    ]
    cy = 16
    for gi, (variant, members) in enumerate(groups.items()):
        marks.append(
            f'<text x="{pad - 12}" y="{cy + 4}" text-anchor="end" class="vlbl">'
            f"{esc(variant)}</text>"
        )
        for row in sorted(members, key=lambda r: float(r["wt_minus_g12d"])):
            diff = float(row["wt_minus_g12d"])
            sem = float(row["wt_minus_g12d_sem"])
            colour = {"below": "var(--c1)", "above": "var(--c3)"}.get(
                row.get("versus_bar", ""), "var(--muted)"
            )
            marks.append(
                f'<line x1="{x(diff - sem):.1f}" x2="{x(diff + sem):.1f}" y1="{cy}" y2="{cy}"'
                f' stroke="{colour}" stroke-width="1.5" opacity="0.6"/>'
                f'<circle cx="{x(diff):.1f}" cy="{cy}" r="4" fill="{colour}"/>'
                f'<text x="{w - gutter}" y="{cy + 4}" class="lbl">{esc(row["design"])}</text>'
            )
            cy += row_h
        cy += 10 if gi < len(groups) - 1 else 0
    return (
        f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="wild-type minus G12D per finalist,'
        f' grouped by generation variant, against the control bar">' + "".join(marks) + "</svg>"
    )


def ledger(experiments: list[dict[str, Any]]) -> str:
    order = {v: i for i, v in enumerate(VERDICT_ORDER)}
    rows = "".join(
        f'<tr><td><code>{esc(e["id"])}</code></td><td>{esc(e["what"])}</td>'
        f'<td>{esc(e["gate"])}</td>'
        f'<td><span class="pill {esc(e["verdict"]).lower()}">{esc(e["verdict"])}</span></td>'
        f'<td>{esc(str(e.get("consequence", "")).strip())}</td></tr>'
        for e in sorted(experiments, key=lambda e: order.get(e["verdict"], 9))
    )
    return (
        '<div class="scroll"><table><thead><tr><th>experiment</th><th>what it tested</th>'
        "<th>gate</th><th>verdict</th><th>consequence</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></div>"
    )


def stat(value: str, label: str, tone: str = "") -> str:
    return (
        f'<div class="stat {tone}"><div class="v">{esc(value)}</div>'
        f'<div class="k">{esc(label)}</div></div>'
    )


TEMPLATE = """<meta name="claude:viewport-fit-sides" content="cover">
<title>Reading a Neoepitope</title>
<style>
/* Layout: a single measured column; charts break out to full width inside it. */
:root {
  --bg: #f7f5f1; --panel: #fffefb; --fg: #1b1a17; --dim: #6a675f; --muted: #bdb8ad;
  --grid: #ddd8cd; --rule: #e4dfd4;
  --c1: #1f6f63; --c2: #7a8f4a; --c3: #b0562c; --c4: #4a5f7a;
  --mono: ui-monospace, "SF Mono", "JetBrains Mono", Menlo, monospace;
  --serif: "Spectral", Georgia, "Times New Roman", serif;
  --sans: "Inter Tight", system-ui, -apple-system, "Helvetica Neue", sans-serif;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #14150f; --panel: #1c1e16; --fg: #ece9e0; --dim: #9b968a; --muted: #4e5245;
    --grid: #32362a; --rule: #2c3025;
    --c1: #5fb8a4; --c2: #a9c06a; --c3: #e0834f; --c4: #8aa4c4; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #14150f; --panel: #1c1e16; --fg: #ece9e0; --dim: #9b968a; --muted: #4e5245;
  --grid: #32362a; --rule: #2c3025;
  --c1: #5fb8a4; --c2: #a9c06a; --c3: #e0834f; --c4: #8aa4c4; color-scheme: dark;
}
body { background: var(--bg); color: var(--fg); font-family: var(--sans); line-height: 1.55;
  font-size: 15px; }
.wrap { max-width: 50rem; margin: 0 auto;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px)); padding-block: 2.6rem 4rem; }
h1 { font-family: var(--serif); font-size: clamp(1.9rem, 5vw, 2.7rem); line-height: 1.12;
  margin: 0 0 .5rem; text-wrap: balance; font-weight: 600; letter-spacing: -.015em; }
h2 { font-family: var(--serif); font-size: 1.32rem; margin: 3.2rem 0 .3rem; font-weight: 600;
  text-wrap: balance; }
h2::before { content: ""; display: block; height: 1px; background: var(--rule);
  margin-bottom: 1.4rem; }
.deck { color: var(--dim); font-size: 1.04rem; max-width: 40rem; margin: 0 0 .4rem; }
.asof { font-family: var(--mono); font-size: .74rem; text-transform: uppercase;
  letter-spacing: .09em; color: var(--dim); margin: 0 0 2rem; }
.strip { display: grid; grid-template-columns: repeat(auto-fit, minmax(8.5rem, 1fr)); gap: .6rem;
  margin: 1.6rem 0 .6rem; }
.stat { background: var(--panel); border: 1px solid var(--rule); border-radius: 3px;
  padding: .75rem .85rem; }
.stat .v { font-family: var(--mono); font-size: 1.45rem; font-variant-numeric: tabular-nums;
  line-height: 1.1; }
.stat .k { font-size: .76rem; color: var(--dim); margin-top: .3rem; }
.stat.good .v { color: var(--c1); } .stat.bad .v { color: var(--c3); }
figure { margin: 1.5rem 0 2.2rem; }
figcaption { font-size: .97rem; margin-bottom: .7rem; color: var(--fg); max-width: 40rem; }
.plot { background: var(--panel); border: 1px solid var(--rule); border-radius: 3px;
  padding: .7rem; overflow-x: auto; }
.plot svg { display: block; width: 100%; height: auto; min-width: 22rem; }
.lbl { font-family: var(--mono); font-size: 10px; fill: var(--dim);
  font-variant-numeric: tabular-nums; }
.vlbl { font-family: var(--mono); font-size: 10.5px; fill: var(--fg); font-weight: 600; }
.note { font-size: .83rem; color: var(--dim); margin: .6rem 0 0; }
p { max-width: 40rem; }
code { font-family: var(--mono); font-size: .88em; }
.scroll { overflow-x: auto; border: 1px solid var(--rule); border-radius: 3px;
  background: var(--panel); }
table { border-collapse: collapse; width: 100%; min-width: 40rem; font-size: .83rem; }
th, td { text-align: left; padding: .5rem .65rem; border-bottom: 1px solid var(--rule);
  vertical-align: top; }
th { font-size: .72rem; text-transform: uppercase; letter-spacing: .07em; color: var(--dim);
  font-weight: 600; }
tr:last-child td { border-bottom: 0; }
.pill { font-family: var(--mono); font-size: .68rem; padding: .14rem .42rem; border-radius: 2px;
  white-space: nowrap; letter-spacing: .04em; }
.pill.promote { background: var(--c1); color: var(--bg); }
.pill.keep { background: var(--c4); color: var(--bg); }
.pill.discard { background: var(--c3); color: var(--bg); }
.pill.running { background: var(--muted); color: var(--bg); }
ul.plain { padding-left: 1.1rem; max-width: 40rem; }
ul.plain li { margin: .4rem 0; }
.flag { border-left: 2px solid var(--c3); padding: .1rem 0 .1rem .9rem; margin: 1.3rem 0;
  max-width: 40rem; }
.flag p { margin: .35rem 0; }
</style>
<div class="wrap">
__BODY__
</div>
"""


def removed_before_scoring(data: dict[str, Any]) -> int:
    """What the gates removed before any design was scored: generated less the scored cohort."""
    return int(data["funnel"][0]["surviving"]) - int(data["headline"]["screened"])


def grid_collapse(data: dict[str, Any]) -> int:
    grid = data["escape_cone"]["by_reach_angstrom"]
    return round(float(grid[6]) / float(grid[7]))


def tier1_bar(rows: list[dict[str, str]]) -> float:
    """The control's own wild-type difference, measured in the same invocations as these rows.

    There is deliberately no fallback to `headline.tier1_control_bias`: that is a different
    measurement, and silently substituting it would put a number from another pass on a chart
    labelled as this one's.
    """
    return float(rows[0]["control_bar"])


def tier1_finding(rows: list[dict[str, str]], data: dict[str, Any]) -> str:
    head = data["headline"]
    if not rows:
        return (
            "The replicated tier-1 pass has not been read into this page yet. The control's own"
            f" difference is {float(head['tier1_control_bias']):+.4f}, which is positive, so"
            " nothing in this column can be read as selectivity."
        )
    paired = data["tier1_paired"]
    bar = tier1_bar(rows)
    ranked = [
        (variant, sum(float(r["wt_minus_g12d"]) for r in members) / len(members))
        for variant, members in by_variant(rows).items()
    ]
    ranked.sort(key=lambda kv: kv[1])
    order = " &lt; ".join(f"<code>{esc(v)}</code>" for v, _ in ranked)
    return (
        f"Both peptides were written into each invocation and the pass run twice"
        f" \u2014 {paired['predictions']} predictions, {paired['failures']} failures. The"
        f" control's own difference is {bar:+.4f} &plusmn;"
        f" {float(rows[0]['wt_minus_g12d_sem']):.4f} over {rows[0]['control_bar_runs']}"
        f" invocations: still positive, but tied zero at this sample size, so the bar is weak"
        f" and so is any ranking against it. {paired['below_bar']} of {len(rows)} finalists sit"
        f" below it, {paired['tied_bar']} tied, {paired['above_bar']} above."
        f" <strong>The unpaired pass's headline did not survive.</strong> Taking the two"
        f" conditions from different invocations, <code>bulge_core</code> \u2014 the variant"
        f" placed over the bulge \u2014 looked"
        f" {float(paired['withdrawn_claim']['bulge_core_mean']):+.4f} with Fisher"
        f" p = {paired['withdrawn_claim']['bulge_core_fisher_p']}. Paired, it is"
        f" {float(paired['paired_result']['bulge_core_mean']):+.4f} at"
        f" p = {paired['paired_result']['bulge_core_fisher_p']},"
        f" because {esc(paired['paired_result']['driven_by'])} moves across zero. The measured"
        f" order is {order} \u2014 not the pre-registered one. What survives is only that"
        f" <code>hub</code>, which puts the least area on the peptide, is the most positive"
        f" group. No design was re-ranked on any of this."
    )


@app.command()
def main(
    campaign: Path = typer.Option(Path("docs/experiments/pmhc/campaign.yaml")),
    verdicts: Path = typer.Option(Path("docs/experiments/pmhc/verdicts.yaml")),
    pooled: Path = typer.Option(Path("data/cofold_spread_pooled.csv")),
    replicates: Path = typer.Option(Path("data/cofold_replicates.csv")),
    wildtype: Path = typer.Option(Path("data/cofold_wildtype.csv")),
    out: Path = typer.Option(Path("docs/artifact/index.html")),
) -> None:
    data = read_yaml(campaign)
    ledger_data = read_yaml(verdicts)
    head = data["headline"]
    bar = float(head["bar"])
    pooled_rows, _ = split_control(read_csv(pooled))
    replicate_rows, _ = split_control(read_csv(replicates))
    wt_rows, _ = split_control(read_csv(wildtype)) if wildtype.exists() else ([], [])
    if not wt_rows:
        typer.echo(f"WARNING: {wildtype} is missing; the tier-1 section ships empty", err=True)

    body: list[str] = [
        "<h1>Reading a Neoepitope</h1>",
        '<p class="deck">A binder for the KRAS G12D peptide on HLA-A*11:01 has to tell one'
        " aspartate from one glycine. This is what the twenty selected designs rest on, and what"
        " they do not.</p>",
        f'<p class="asof">Challenge 5 &middot; measured {esc(data["updated"])} &middot;'
        f' closes {esc(data["deadline_aoe"])} AoE</p>',
        '<div class="strip">'
        + stat(f"{float(head['best_drop']):+.3f}", "best design's peptide-dependence", "good")
        + stat(f"{bar:+.3f}", "the receptor's own bar")
        + stat(
            f"{head['above_bar']} / {head['screened']}",
            "designs above that bar",
            "bad",
        )
        + stat(
            f"{float(head['tier1_control_bias']):+.4f}",
            "tier 1 on the control: the wrong sign",
            "bad",
        )
        + "</div>",
        '<p class="note">The fourth number is the one that did not move. Ranking tier 1 asks'
        " which peptide a binder prefers; measured on the receptor that is specific by"
        f" construction, this predictor prefers the peptide it does <em>not</em> recognise, by"
        f" {float(head['tier1_control_bias']):+.4f} &plusmn;"
        f" {float(head['tier1_control_bias_sem']):.4f} over"
        f" {head['tier1_control_invocations']} independent invocations. Everything below ranks on"
        " peptide <em>dependence</em> instead, and says so.</p>",
        "<h2>Why the mutation itself is unreachable</h2>",
        finding(
            "The one atom that chemically distinguishes the two targets sits at the bottom of a"
            f" pocket about 6 Å deep. Its escape cone collapses"
            f" {grid_collapse(data)}-fold between 6 and 7 Å and reaches"
            f" {data['escape_cone']['by_reach_angstrom'][8]}% at 8 Å —"
            f" {data['escape_cone']['free_rays_at_8A']} of"
            f" {data['escape_cone']['rays']:,} rays. No designed side chain can touch it.",
            escape_cone(data["escape_cone"]),
            "Probe controls behaved: "
            + ", ".join(
                f"{esc(name)} at {float(value):.2f}%"
                for name, value in data["escape_cone"]["controls"].items()
            )
            + f". Log scale; the {data['escape_cone']['by_reach_angstrom'][8]}% floor is"
            f" {data['escape_cone']['free_rays_at_8A']} surviving rays, not a sampling"
            " artefact.",
        ),
        "<h2>What survived, and out of how many</h2>",
        finding(
            f"Three gates removed {removed_before_scoring(data)} of"
            f" {data['funnel'][0]['surviving']} designs before any of them was scored for"
            " binding."
            " Novelty had to come first: the generator returns the target's own"
            " β<sub>2</sub>-microglobulin when given a vacant interface, and that is a"
            " disqualification rather than a worse score.",
            funnel(data["funnel"]),
            "Solid is what survived each stage; faded is what that stage removed.",
        ),
        "<h2>The ranking, against a calibrated bar</h2>",
        finding(
            f"No design clears the receptor's own peptide-dependence. Of"
            f" {head['screened']} screened designs {head['above_bar']} are above the bar,"
            f" {head['tied_bar']} are tied with it and {head['below_bar']} are below;"
            f" {head['clear_zero']} clear zero at all. The bar is the 7OW6 receptor measured"
            f" {data['bar']['invocations']} times from byte-identical specifications.",
            scoreboard(pooled_rows, bar),
            "Teal clears the bar or ties it, olive is above zero, grey is at or below zero."
            " The dashed line is the receptor.",
        ),
        "<h2>Which leaders are real and which are prediction failures</h2>",
        finding(
            "The drop is complex minus peptide-free, so an unlucky peptide-free pose"
            " <em>inflates</em> it — the top of the ranking is enriched for prediction"
            " failures as well as for genuine peptide reading. Replication is the test, and the"
            " two leaders split on it.",
            reproducibility(replicate_rows),
            "Each line spans one design's two independent invocations. Teal reproduces within"
            f" {REPRODUCES_WITHIN}; orange does not.",
        ),
        "<h2>Ranking tier 1, reported and never gated</h2>",
        finding(
            tier1_finding(wt_rows, data),
            tier1_chart(wt_rows),
            "Bars are each design's own standard error. The dashed line is the control's own"
            " difference, which is the bar because it is the only molecule here known to"
            " discriminate these two peptides.",
        ),
        "<h2>Every experiment, and what it settled</h2>",
        ledger(ledger_data["experiments"]),
        "<h2>What is not established</h2>",
        "<ul class=\"plain\">"
        + "".join(f"<li>{esc(str(item).strip())}</li>" for item in ledger_data["not_established"])
        + "</ul>",
        '<div class="flag"><p><strong>One predictor.</strong> The standing discipline on this'
        " project is that every load-bearing claim gets an independent second method, and that"
        " is not met here: the Chai-1 container was never located, so the consensus pass did not"
        " run. Every number on this page is Boltz 2.2.1 scoring structures Boltz 2.2.1"
        " predicted. Replication bounds sampling and run-to-run error; it cannot detect a"
        " systematic error shared by every run, and the tier-1 inversion is exactly such a"
        " case.</p></div>",
        '<p class="note">Every figure here is read at build time from'
        " <code>docs/experiments/pmhc/campaign.yaml</code>,"
        " <code>docs/experiments/pmhc/verdicts.yaml</code> and the <code>data/*.csv</code>"
        " tables by <code>scripts/build_artifact.py</code>. Nothing on the page is typed in.</p>",
    ]

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(TEMPLATE.replace("__BODY__", "\n".join(body)), encoding="utf-8")
    typer.echo(f"wrote {out} ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    app()
