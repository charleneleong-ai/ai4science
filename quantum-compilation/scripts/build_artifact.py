"""Render the VQCC progress artifact from the verifiers' own outputs.

Numbers come from three committed files -- the official verifier runs over
`sample_submission/` in `tests/data/reference_<task>.json`, the verdict ledger in
`docs/experiments/verdicts.yaml`, and the dated public-leaderboard snapshot in
`docs/experiments/leaderboard-<date>.json` -- plus the formulas and constants imported from
`quantum_compilation.cost`, so the page cannot drift from the module the project scores
with. `HAND_RECORDED` holds the handful of figures that exist in no file (wall-clock run
times, the competition's submission budget); the footer says so.

    uv run --extra tracking python scripts/build_artifact.py

Writes `docs/artifact/index.html`, committed on every rebuild and republished as a Claude
Artifact at each significant result update.
"""

from __future__ import annotations

import html
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from quantum_compilation.cost import (
    MAX_ANGLE_QUBITS,
    MAX_DQI_FAILURE_RATE,
    MAX_HAMILTONIAN_ERROR_HA,
    MIN_E_GAP_HA,
    TASKS,
    MetricCalibration,
    dqi_repetition_factor,
    rotation_t_count,
    score,
    task_ratio,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = PROJECT_ROOT / "tests" / "data"
EXPERIMENTS_DIR = PROJECT_ROOT / "docs" / "experiments"
VERDICTS_PATH = EXPERIMENTS_DIR / "verdicts.yaml"
LEADERBOARD_PATH = EXPERIMENTS_DIR / "leaderboard-2026-10-31.json"
OUTPUT_PATH = PROJECT_ROOT / "docs" / "artifact" / "index.html"

STATUS_DATE = "2026-10-31"
SUBMISSION_ID = "57824925"
OUR_SCORE = 1.00000

# Wall-clock minutes per verifier on 10 CPU cores, and the competition's daily submission
# budget. Neither is in any result file; the footer declares them as recorded by hand.
HAND_RECORDED = {"minutes": {"invsqrt": 5.7, "dqi": 2.4, "dfthcsos": 1.0}, "submissions_per_day": 5}

# Self-reported by a competitor in Kaggle discussion topic 754174 on 2026-10-30 and
# re-verified by them, not by us. Five significant figures, which bounds how finely the
# blend cross-check below can agree. Labelled third-party wherever it appears.
PUBLIC_COSTS = {"invsqrt": 100_941.0, "dqi": 4.0334e11, "dfthcsos": 3.3780e8}
PUBLIC_BLEND = 23.13587
# Half of the last quoted unit in each figure above, as written -- 100941 to the unit,
# 4.0334e11 and 3.3780e8 to five significant figures. This bounds how finely the blend
# cross-check can agree, so it is stated rather than inferred from the float.
PUBLIC_COST_HALF_ULP = {"invsqrt": 0.5, "dqi": 5e6, "dfthcsos": 5e3}
# Shown exactly as the competitor quoted them, so the table cannot imply digits they did not
# report (a plain thousands format turns 3.3780e8 into a 9-figure 337,800,000).
PUBLIC_COST_QUOTED = {"invsqrt": "100,941", "dqi": "4.0334×10¹¹", "dfthcsos": "3.3780×10⁸"}

# How each task's cost is assembled, and which reported field carries its leading factor.
# `cost = factor * T^t_power * Q^q_power`, verified against cost.py's three cost functions.
SCORED_AS = {
    "invsqrt": ("Q × C", "toffoli_equivalent", "C", 0, 1),
    "dqi": ("R(ε₊) × Q² × T", "repetition_factor", "R", 1, 2),
    "dfthcsos": ("λ_eff × T × Q", "lambda_eff", "λ_eff", 1, 1),
}

SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
INLINE_CODE = re.compile(r"`([^`]+)`")
INLINE_BOLD = re.compile(r"\*\*([^*]+)\*\*")


@dataclass(frozen=True)
class TaskReport:
    """One task's verifier output, with the handful of fields the page draws."""

    task: str
    detail: dict[str, Any]

    @property
    def cost(self) -> float:
        return float(self.detail["score"])

    @property
    def qubits(self) -> int:
        return int(self.detail["num_qubits"])

    @property
    def toffolis(self) -> float:
        return float(self.detail["toffoli_count"])

    @property
    def shots(self) -> int:
        return int(self.detail["verified_shots"])

    @property
    def eps_plus(self) -> float:
        return float(self.detail["effective_error_rate"])

    @property
    def repetition_factor(self) -> float:
        return float(self.detail["repetition_factor"])

    @property
    def leading_factor(self) -> tuple[str, str, float]:
        """The cost's written form, its leading factor's name, and that factor's value."""
        formula, field_name, symbol, _, _ = SCORED_AS[self.task]
        return formula, symbol, float(self.detail[field_name])


def load_reports() -> dict[str, TaskReport]:
    return {
        task: TaskReport(task, json.loads((REFERENCE_DIR / f"reference_{task}.json").read_text()))
        for task in TASKS
    }


def inline_markup(text: str) -> str:
    """Escape YAML prose, then honour the `code` and **bold** the ledger is written in."""
    escaped = html.escape(text.strip(), quote=False)
    escaped = INLINE_CODE.sub(r"<code>\1</code>", escaped)
    return INLINE_BOLD.sub(r"<strong>\1</strong>", escaped)


def fmt(value: float) -> str:
    """A number a reader can compare at a glance, with its magnitude kept visible."""
    if value < 1e10:
        return f"{value:,.2f}"
    exponent = int(math.floor(math.log10(value)))
    mantissa = value / 10**exponent
    if mantissa >= 9.999995:  # would render as 10.00000x10^n
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"{mantissa:.5f}×10{str(exponent).translate(SUPERSCRIPT)}"


def blend_precision_ppm(costs: dict[str, float], half_ulp: dict[str, float]) -> float:
    """How finely a blend of rounded costs can agree, in parts per million.

    Each cost's quoted rounding is a relative error; a geometric mean over n tasks scales
    each by 1/n, and they add in quadrature.
    """
    relative = [half_ulp[task] / cost for task, cost in costs.items()]
    return math.sqrt(sum((r / len(costs)) ** 2 for r in relative)) * 1e6


def svg(width: int, height: int, aria: str, body: str) -> str:
    return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{aria}">{body}</svg>'


def hgrid(x0: float, x1: float, y: float, label: str) -> str:
    return (
        f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" class="grid"/>'
        f'<text x="{x0 - 10:.1f}" y="{y + 6:.1f}" text-anchor="end" class="ct">{label}</text>'
    )


def field_svg(scores: list[float], rank_best: int, rank_worst: int) -> str:
    """The whole field, sorted: the plateau, the cliff, and which block we are in."""
    width, height = 680, 400
    pad_l, pad_r, pad_t, pad_b = 74, 26, 24, 62
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b
    teams = len(scores)
    y_hi = max(scores) * 1.06

    def place(rank: int, value: float) -> tuple[float, float]:
        return (
            pad_l + (rank - 1) / (teams - 1) * plot_w,
            pad_t + plot_h - value / y_hi * plot_h,
        )

    path = []
    for rank, value in enumerate(scores, start=1):
        x, y = place(rank, value)
        path.append(f"{'M' if not path else 'L'}{x:.1f} {y:.1f}")

    grid = "".join(
        hgrid(pad_l, pad_l + plot_w, place(1, tick)[1], str(tick))
        for tick in (0, 5, 10, 15, 20, 25)
    )
    ticks = "".join(
        f'<text x="{place(rank, 0)[0]:.1f}" y="{height - 34}" text-anchor="middle" '
        f'class="ct">{rank:,}</text>'
        for rank in (1, 500, 1000, 1500, teams)
    )
    x_best, x_worst = place(rank_best, 0)[0], place(rank_worst, 0)[0]
    our_y = place(rank_best, OUR_SCORE)[1]
    leader_x, leader_y = place(1, scores[0])
    return svg(
        width,
        height,
        f"All {teams:,} public leaderboard scores in rank order: a dense plateau near 23, a "
        f"cliff, and our score of 1.0 among the lowest",
        grid
        + f'<rect x="{x_best:.1f}" y="{pad_t}" width="{max(3.0, x_worst - x_best):.1f}" '
        f'height="{plot_h}" fill="var(--warn)" fill-opacity="0.18"/>'
        + ticks
        + f'<path d="{" ".join(path)}" fill="none" stroke="var(--accent)" stroke-width="3"/>'
        + f'<circle cx="{leader_x:.1f}" cy="{leader_y:.1f}" r="7" fill="var(--accent)"/>'
        + f'<text x="{leader_x + 14:.1f}" y="{leader_y + 7:.1f}" class="cv">'
        f"leader {scores[0]:.2f}</text>"
        + f'<text x="{x_worst - 10:.1f}" y="{our_y - 18:.1f}" text-anchor="end" class="cv">'
        f"us {OUR_SCORE:.2f}</text>"
        + f'<text x="{x_worst - 10:.1f}" y="{our_y + 2:.1f}" text-anchor="end" class="ct">'
        f"ranks {rank_best:,}–{rank_worst:,}</text>"
        + f'<text x="{pad_l + plot_w / 2:.0f}" y="{height - 8}" text-anchor="middle" '
        f'class="cax">rank, of {teams:,} teams</text>'
        + f'<text x="22" y="{pad_t + plot_h / 2:.0f}" class="cax" '
        f'transform="rotate(-90 22 {pad_t + plot_h / 2:.0f})" text-anchor="middle">score</text>',
    )


def repetition_svg(dqi: TaskReport) -> str:
    """R(eps) against the measured bound: the 19% the reference pays for zero failures."""
    width, height = 680, 330
    pad_l, pad_r, pad_t, pad_b = 66, 30, 26, 62
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b
    y_hi = 3.0
    x_lo, x_hi = math.log10(1e-5), math.log10(MAX_DQI_FAILURE_RATE)

    def place(eps: float, r: float) -> tuple[float, float]:
        return (
            pad_l + (math.log10(eps) - x_lo) / (x_hi - x_lo) * plot_w,
            pad_t + plot_h - min(r, y_hi) / y_hi * plot_h,
        )

    path = []
    for step in range(241):
        eps = 10 ** (x_lo + (x_hi - x_lo) * step / 240)
        rate = dqi_repetition_factor(eps)
        if rate > y_hi:
            break
        x, y = place(eps, rate)
        path.append(f"{'M' if not path else 'L'}{x:.1f} {y:.1f}")

    grid = "".join(
        hgrid(pad_l, pad_l + plot_w, place(1e-5, tick)[1], f"{tick}×") for tick in (1, 2, 3)
    )
    ticks = "".join(
        f'<text x="{place(10.0**decade, 0)[0]:.1f}" y="{height - 34}" text-anchor="middle" '
        f'class="ct">10{str(decade).translate(SUPERSCRIPT)}</text>'
        for decade in (-5, -4, -3, -2)
    )
    mx, my = place(dqi.eps_plus, dqi.repetition_factor)
    pole_x = pad_l + plot_w
    return svg(
        width,
        height,
        "The dqi repetition factor against the bounded failure rate, with the reference's "
        "1.19 times marked and the rejection pole at one sixteenth",
        grid
        + ticks
        + f'<line x1="{pole_x:.1f}" y1="{pad_t}" x2="{pole_x:.1f}" y2="{pad_t + plot_h}" '
        f'class="pole"/>'
        + f'<text x="{pole_x - 10:.1f}" y="{pad_t + 20}" text-anchor="end" class="cl2">'
        f"rejected at ε = 1/16</text>"
        + f'<path d="{" ".join(path)}" fill="none" stroke="var(--accent)" stroke-width="3"/>'
        + f'<line x1="{pad_l}" y1="{my:.1f}" x2="{mx:.1f}" y2="{my:.1f}" class="marker"/>'
        + f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="8" fill="var(--warn)"/>'
        + f'<text x="{mx + 14:.1f}" y="{my - 4:.1f}" class="cv">'
        f"R = {dqi.repetition_factor:.4f}</text>"
        + f'<text x="{pad_l + plot_w / 2:.0f}" y="{height - 8}" text-anchor="middle" '
        f'class="cax">bounded failure rate ε₊ (log)</text>',
    )


def angle_share_svg(report: TaskReport) -> str:
    """Two true partitions: what the ANGLE circuits are of the Toffolis, and of the width.

    Only additive decompositions go in a part-of-whole bar. The scored total width (730) is
    *not* one -- `temporary_qubits` is a max over the three families, not a sum -- so the
    bars show the Toffoli partition and the persistent-qubit partition, and the prose below
    carries the max.
    """
    resources = report.detail["resources"]
    bars = [
        (
            "Toffolis",
            float(resources["angles_toffoli"]),
            float(resources["toffoli"]),
            "of all Toffolis",
        ),
        (
            "Persistent qubits",
            float(resources["angles_persistent_qubits"]),
            float(resources["persistent_qubits"]),
            "of persistent width",
        ),
    ]
    width, bar_h, gap = 680, 52, 46
    pad_l, pad_r = 196, 200
    plot = width - pad_l - pad_r
    height = len(bars) * (bar_h + gap) + 22
    rows = []
    for index, (label, part, whole, caption) in enumerate(bars):
        y = index * (bar_h + gap) + 18
        frac = part / whole
        rows.append(
            f'<text x="{pad_l - 16}" y="{y + 32}" text-anchor="end" class="cl">{label}</text>'
            f'<rect x="{pad_l}" y="{y}" width="{plot}" height="{bar_h}" rx="5" '
            f'fill="var(--neutral-mark)" fill-opacity="0.3"/>'
            f'<rect x="{pad_l}" y="{y}" width="{frac * plot:.1f}" height="{bar_h}" rx="5" '
            f'fill="var(--accent)"/>'
            f'<text x="{pad_l + 14}" y="{y + 33}" class="cvi">{frac * 100:.0f}%</text>'
            f'<text x="{pad_l + plot + 14}" y="{y + 25}" class="cv">'
            f"{part:,.0f} / {whole:,.0f}</text>"
            f'<text x="{pad_l + plot + 14}" y="{y + 43}" class="ct">{caption}</text>'
        )
    return svg(
        width,
        height,
        "The ANGLE circuits are most of the chemistry task's Toffoli count but about half of "
        "its persistent qubits",
        "".join(rows),
    )


def reference_rows(reports: dict[str, TaskReport]) -> str:
    rows = []
    for task in TASKS:
        report = reports[task]
        formula, symbol, value = report.leading_factor
        rows.append(
            f"<tr><td><code>{task}</code></td>"
            f'<td class="num">{fmt(report.cost)}</td>'
            f"<td><code>{formula}</code></td>"
            f'<td class="num">{report.qubits:,}</td>'
            f'<td class="num">{report.toffolis:,.1f}</td>'
            f'<td class="num">{symbol} = {value:,.4f}</td>'
            f'<td class="num">{report.shots:,}</td></tr>'
        )
    return "".join(rows)


def ratio_rows(reports: dict[str, TaskReport], calibration: MetricCalibration) -> str:
    return "".join(
        f"<tr><td><code>{task}</code></td>"
        f'<td class="num">{PUBLIC_COST_QUOTED[task]}</td>'
        f'<td class="num">{task_ratio(task, PUBLIC_COSTS[task], calibration):,.3f}×</td></tr>'
        for task in TASKS
    )


def run_log_rows(reports: dict[str, TaskReport]) -> str:
    return "".join(
        f'<div class="fcell"><div class="ftask">{task}</div>'
        f'<div class="fpass">PASS</div>'
        f'<div class="fmeta">{reports[task].shots:,} shots</div>'
        f'<div class="fmeta">{HAND_RECORDED["minutes"][task]} min · '
        f"cost {fmt(reports[task].cost)}</div></div>"
        for task in TASKS
    )


def verdict_blocks(verdicts: dict[str, Any]) -> str:
    tone = {"PROMOTE": "ok", "KEEP": "warn", "DISCARD": "bad", "RUNNING": "run"}
    blocks = []
    for item in verdicts.get("experiments", []):
        verdict = str(item.get("verdict", "RUNNING"))
        pr = item.get("pr")
        link = (
            f' · <a href="https://github.com/charleneleong-ai/ai4science/pull/{pr}">PR #{pr}</a>'
            if pr
            else ""
        )
        blocks.append(
            f'<article class="verdict"><header>'
            f"<code>{html.escape(str(item['name']), quote=False)}</code>"
            f'<span class="pill {tone.get(verdict, "run")}">{html.escape(verdict)}</span>'
            f'</header><p class="vmeta">gate {html.escape(str(item.get("gate", "—")))}'
            f" · {html.escape(str(item.get('date', '')))}{link}</p>"
            f"<p>{inline_markup(str(item.get('result', '')))}</p></article>"
        )
    return "".join(blocks)


def gate_rows(verdicts: dict[str, Any]) -> str:
    return "".join(
        f"<tr><td><code>{html.escape(name)}</code></td><td>{inline_markup(text)}</td></tr>"
        for name, text in verdicts.get("gates", {}).items()
    )


def build() -> str:
    reports = load_reports()
    verdicts = yaml.safe_load(VERDICTS_PATH.read_text())
    board = json.loads(LEADERBOARD_PATH.read_text())
    scores: list[float] = board["scores"]
    teams = board["teams"]

    invsqrt, dqi, dfthcsos = (reports[t] for t in TASKS)
    calibration = MetricCalibration(
        reference_costs={t: reports[t].cost for t in TASKS},
        provenance=f"measured by the shipped verifiers over sample_submission/, {STATUS_DATE}",
        measured=True,
    )

    above = sum(1 for s in scores if s > OUR_SCORE + 5e-6)
    tied = sum(1 for s in scores if abs(s - OUR_SCORE) < 5e-6)
    rank_best, rank_worst = above + 1, above + tied
    median = scores[teams // 2]
    plateau_hi = scores[int(teams * 0.10)]
    escaped = sum(1 for s in scores if s >= 25.0)
    leader = scores[0]
    gap_score = leader / OUR_SCORE
    gap_cost = gap_score**3

    blend = score(PUBLIC_COSTS, calibration)
    blend_delta_ppm = abs(blend - PUBLIC_BLEND) / PUBLIC_BLEND * 1e6
    blend_floor_ppm = blend_precision_ppm(PUBLIC_COSTS, PUBLIC_COST_HALF_ULP)

    rotation_count = int(invsqrt.detail["rotation_count"])
    rotation_eps = float(invsqrt.detail["rotation_eps"])
    t_per_rotation = rotation_t_count(rotation_eps)
    resources = dfthcsos.detail["resources"]
    metrics = dfthcsos.detail["metrics"]
    angle_qubits = int(resources["angles_qubits"])
    error_2 = float(metrics["error_2"])
    e_gap = float(metrics["E_gap"])
    toffoli_span = dqi.toffolis / invsqrt.toffolis

    return f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>VQCC Compression Log</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>
/* Instrument readout: one measured column, data in mono, the finding above each chart. */
:root {{
  --bg: #f7f8f7; --surface: #ffffff; --fg: #15201c; --muted: #5b6b64;
  --line: #dfe4e1; --grid: #b4c0ba; --accent: #0f6f5c; --accent-soft: #e4f0ec;
  --warn: #b5700f; --warn-soft: #fbeedd; --bad: #a3342a; --bad-soft: #f8e5e2;
  --neutral-mark: #8a9890;
  --display: "IBM Plex Serif", Georgia, serif;
  --body: "IBM Plex Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #0e1412; --surface: #151d1a; --fg: #e6ece9; --muted: #93a39c;
    --line: #26322e; --grid: #46554f; --accent: #4fc4a6; --accent-soft: #16302a;
    --warn: #e0a449; --warn-soft: #33270f; --bad: #e07a6e; --bad-soft: #351e1b;
    --neutral-mark: #5d6d67;
    color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --bg: #0e1412; --surface: #151d1a; --fg: #e6ece9; --muted: #93a39c;
  --line: #26322e; --grid: #46554f; --accent: #4fc4a6; --accent-soft: #16302a;
  --warn: #e0a449; --warn-soft: #33270f; --bad: #e07a6e; --bad-soft: #351e1b;
  --neutral-mark: #5d6d67;
  color-scheme: dark;
}}
* {{ box-sizing: border-box; }}
body {{
  background: var(--bg); color: var(--fg); font-family: var(--body);
  font-size: 15px; line-height: 1.6; margin: 0;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 40px 72px;
}}
main {{ max-width: 760px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }}
h1 {{ font-family: var(--display); font-size: clamp(27px, 5vw, 38px); line-height: 1.15;
      margin: 0 0 12px; text-wrap: balance; font-weight: 600; }}
h2 {{ font-family: var(--display); font-size: 21px; margin: 0 0 4px; font-weight: 600; }}
h3 {{ font-size: 13px; text-transform: uppercase; letter-spacing: 0.08em;
      color: var(--muted); margin: 0 0 10px; font-weight: 600; }}
p {{ margin: 0 0 12px; }}
section {{ display: flex; flex-direction: column; gap: 14px; }}
.status {{ font-family: var(--mono); font-size: 12.5px; color: var(--muted); }}
.lede {{ font-size: 18px; line-height: 1.5; }}
.strip {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; }}
.tile {{ background: var(--surface); border: 1px solid var(--line); border-radius: 8px;
         padding: 14px 16px; display: flex; flex-direction: column; gap: 3px; min-width: 0; }}
.tile .k {{ font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.07em;
            color: var(--muted); font-weight: 600; }}
.tile .v {{ font-family: var(--mono); font-size: 23px; font-weight: 600;
            font-variant-numeric: tabular-nums; }}
.tile .n {{ font-size: 12.5px; color: var(--muted); }}
.tile.flag {{ border-color: var(--warn); background: var(--warn-soft); }}
.finding {{ font-size: 17px; line-height: 1.45; font-weight: 500; margin: 0;
            border-left: 3px solid var(--accent); padding: 2px 0 2px 14px; }}
figure {{ margin: 0; background: var(--surface); border: 1px solid var(--line);
          border-radius: 8px; padding: 16px; display: flex; flex-direction: column; gap: 10px; }}
figure svg {{ width: 100%; height: auto; max-width: 100%; display: block; }}
figcaption {{ font-size: 12.5px; color: var(--muted); }}
.grid {{ stroke: var(--grid); stroke-width: 1; }}
.pole {{ stroke: var(--bad); stroke-width: 2; stroke-dasharray: 5 4; }}
.marker {{ stroke: var(--warn); stroke-width: 1.5; stroke-dasharray: 4 4; }}
text {{ font-family: var(--body); fill: var(--fg); }}
.cl {{ font-size: 20px; fill: var(--fg); }}
.cl2 {{ font-size: 19px; fill: var(--muted); font-weight: 600; }}
.ct {{ font-size: 18px; fill: var(--muted); font-family: var(--mono); }}
.cax {{ font-size: 19px; fill: var(--muted); }}
.cv {{ font-size: 20px; fill: var(--fg); font-family: var(--mono); font-weight: 600; }}
.cvi {{ font-size: 20px; fill: var(--bg); font-family: var(--mono); font-weight: 600; }}
.tw {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 8px;
       background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 13.5px; }}
th, td {{ text-align: left; padding: 9px 12px; border-bottom: 1px solid var(--line);
          vertical-align: top; }}
th {{ font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.06em;
      color: var(--muted); font-weight: 600; white-space: nowrap; }}
tr:last-child td {{ border-bottom: none; }}
td.num {{ font-family: var(--mono); font-variant-numeric: tabular-nums; white-space: nowrap; }}
code {{ font-family: var(--mono); font-size: 0.92em; }}
.pill {{ display: inline-block; font-family: var(--mono); font-size: 11px; font-weight: 600;
         padding: 2px 7px; border-radius: 4px; letter-spacing: 0.04em; }}
.pill.ok {{ background: var(--accent-soft); color: var(--accent); }}
.pill.warn {{ background: var(--warn-soft); color: var(--warn); }}
.pill.bad {{ background: var(--bad-soft); color: var(--bad); }}
.pill.run {{ background: var(--surface); color: var(--muted); border: 1px solid var(--line); }}
.funnel {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; }}
.fcell {{ border: 1px solid var(--line); border-left: 3px solid var(--neutral-mark);
          border-radius: 6px; padding: 12px 14px; background: var(--surface);
          display: flex; flex-direction: column; gap: 2px; min-width: 0; }}
.ftask {{ font-family: var(--mono); font-weight: 600; font-size: 14px; }}
.fpass {{ font-family: var(--mono); font-size: 11px; font-weight: 600; color: var(--muted);
          letter-spacing: 0.05em; }}
.fmeta {{ font-size: 12px; color: var(--muted); font-family: var(--mono); }}
.verdict {{ border: 1px solid var(--line); border-radius: 8px; background: var(--surface);
            padding: 14px 16px; }}
.verdict header {{ display: flex; gap: 10px; align-items: baseline; flex-wrap: wrap; }}
.verdict header code {{ font-size: 13.5px; font-weight: 600; }}
.verdict p {{ margin: 0; font-size: 14px; }}
.verdict .vmeta {{ font-size: 12px; color: var(--muted); font-family: var(--mono);
                   margin: 2px 0 8px; }}
ul {{ margin: 0; padding-left: 20px; }}
li {{ margin-bottom: 8px; }}
li::marker {{ color: var(--neutral-mark); }}
.note {{ border-left: 3px solid var(--neutral-mark); padding: 2px 0 2px 14px;
         color: var(--muted); font-size: 14px; }}
footer {{ border-top: 1px solid var(--line); padding-top: 16px; font-size: 12.5px;
          color: var(--muted); }}
a {{ color: var(--accent); }}
</style>

<main>
  <header>
    <p class="status">GOOGLE QUANTUM AI · VERIFIABLE QUANTUM CIRCUIT COMPILATION · {STATUS_DATE}</p>
    <h1>The metric, measured</h1>
    <p class="lede">The cost model is now read off the three shipped verifiers rather than the
    Evaluation page, and every formula reproduces the verifier&rsquo;s own score exactly. The
    unmodified reference submission scores <strong>1.00000</strong> locally and
    <strong>1.00000</strong> on the leaderboard, confirming the two agree on a known input;
    because the verifiers themselves ship, candidates can then be ranked locally without
    spending submissions. No compression has been attempted, so we sit near the bottom of the
    field.</p>
  </header>

  <section>
    <div class="strip">
      <div class="tile"><span class="k">Our score</span>
        <span class="v">{OUR_SCORE:.5f}</span>
        <span class="n">local = leaderboard · submission {SUBMISSION_ID}</span></div>
      <div class="tile flag"><span class="k">Rank</span>
        <span class="v">{rank_best:,}</span>
        <span class="n">of {teams:,} · tied with {tied - 1} others on the reference</span></div>
      <div class="tile flag"><span class="k">Gap to leader</span>
        <span class="v">×{gap_score:.1f}</span>
        <span class="n">≈{gap_cost:,.0f}× on the product of the three cost ratios</span></div>
      <div class="tile"><span class="k">Compression attempted</span>
        <span class="v">none</span><span class="n">everything here is apparatus</span></div>
    </div>
    <p class="note">1.00000 is the reference&rsquo;s score by definition, not an achievement.
    What this milestone establishes is that the instrument reads true: the thing we will measure
    thousands of candidates with agrees with the leaderboard on a known input.</p>
  </section>

  <section>
    <h2>Where we stand</h2>
    <p class="finding">The field is bimodal. Half of it is packed between {median:.2f} and
    {plateau_hi:.2f} — a band 0.6 wide, almost certainly one shared published approach — and only
    {escaped} teams have escaped past 25. We are in the {tied}-team block at exactly 1.00000 that
    submitted the reference unchanged, ranked {rank_best:,}–{rank_worst:,} of {teams:,}.</p>
    <figure>
      {field_svg(scores, rank_best, rank_worst)}
      <figcaption>Every public score on {STATUS_DATE}, in rank order, from the competition&rsquo;s
      own leaderboard snapshotted in
      <code>docs/experiments/leaderboard-{STATUS_DATE}.json</code>. The shaded band is our tie
      block.</figcaption>
    </figure>
    <p>The crowded plateau is the useful part: a field that dense at one value is a field running
    one parameterisation, and the {escaped} teams above 25 found something it does not reach. That
    is the gap worth understanding.</p>
  </section>

  <section>
    <h2>What the Evaluation page omits</h2>
    <p class="finding">Three constants the page does not state. None changes the cost much — the
    back-derived estimates they correct were within 0.04% — but each changes <em>how</em> cost has
    to be computed, and together they rule out computing it from a circuit at all.</p>
    <div class="tw"><table>
      <thead><tr><th>Quantity</th><th>What the verifier does</th><th>On the reference</th></tr></thead>
      <tbody>
        <tr><td>Rotation price</td>
          <td><code>0.56·log₂(1/ε) + 5.3</code> T gates at one <em>circuit-wide</em> ε, picked
          from the accuracy budget the simulated error leaves over — not a per-Clifford-level
          cost</td>
          <td class="num">{rotation_count} rotations at
          ε={rotation_eps:.2e}<br>{t_per_rotation:,.2f} T gates each</td></tr>
        <tr><td>Toffoli count</td>
          <td>the mean executed <em>per shot</em>, not the static gate count; measurement-based
          uncomputation makes it data-dependent</td>
          <td class="num">dqi executes {dqi.toffolis:,.1f}<br>of
          {int(dqi.detail["max_toffoli_count"]):,} possible</td></tr>
        <tr><td>Failure rate ε₊</td>
          <td>a one-sided Clopper–Pearson <em>upper bound</em> at δ=2⁻¹²⁸, not the observed
          fraction</td>
          <td class="num">0 of {dqi.shots:,} failures<br>bounded at {dqi.eps_plus:.3e}</td></tr>
      </tbody>
    </table></div>
    <p>The consequence is structural: <strong>no task&rsquo;s cost is computable from a static
    reading of a circuit</strong>, <code>invsqrt</code> included. The earlier static pricing path
    was removed rather than corrected — it would have ranked candidates on a number the
    leaderboard does not use.</p>
    <p class="note">The page uses &ldquo;T gates&rdquo; for the rotation price and
    &ldquo;Toffolis&rdquo; for the gate count the cost formulas call T. They are different
    quantities: <code>invsqrt</code> is scored on C = Toffolis + T&nbsp;gates/4.</p>
  </section>

  <section>
    <h2>The 19% nobody can optimise away</h2>
    <p class="finding">A decoder that fails zero of {dqi.shots:,} shots is still bounded at
    ε₊ = {dqi.eps_plus:.2e}, so <code>dqi</code> pays R = {dqi.repetition_factor:.4f} — a 19%
    repetition overhead — before any circuit work. Only more shots lower that floor, and the shot
    count is the verifier&rsquo;s to choose.</p>
    <figure>
      {repetition_svg(dqi)}
      <figcaption>R(ε) = ln2 / −ln(½ + 2√ε), drawn from
      <code>cost.dqi_repetition_factor</code>. The cost diverges at ε = 1/16, where the verifier
      rejects outright. Our point is the Clopper–Pearson bound for 0 failures in
      {dqi.shots:,} shots.</figcaption>
    </figure>
  </section>

  <section>
    <h2>What a chemistry submission actually controls</h2>
    <p class="finding">The ANGLE circuits a submission ships are 85% of the Toffoli count but only
    about half the persistent width — and they add nothing at all to the temporary workspace,
    which is a <em>maximum</em> over families and is dominated by inner state prep. Width is
    governed by a parameterisation, not by the circuits a submission writes.</p>
    <figure>
      {angle_share_svg(dfthcsos)}
      <figcaption>Both bars are additive partitions:
      {resources["angles_toffoli"]:,.0f} + {resources["heuristic_toffoli"]:,} Toffolis, and
      {resources["angles_persistent_qubits"]} + {resources["outer_state_prep_persistent_qubits"]}
      + {resources["inner_state_prep_persistent_qubits"]} persistent qubits. The scored total
      width of {int(resources["qubits"])} is <em>not</em> a partition —
      {int(resources["persistent_qubits"])} persistent +
      {int(resources["temporary_qubits"])} temporary, where the temporary term is
      max({resources["angles_temporary_qubits"]},
      {resources["outer_state_prep_temporary_qubits"]},
      {resources["inner_state_prep_temporary_qubits"]}).</figcaption>
    </figure>
    <p class="note">The reference sits close to three cliffs, all of which reject outright rather
    than gracefully: the widest ANGLE circuit at {angle_qubits} qubits against the
    {MAX_ANGLE_QUBITS} cap ({angle_qubits / MAX_ANGLE_QUBITS * 100:.0f}% of it),
    <code>error_2</code> at {error_2:.5f} Ha against a {MAX_HAMILTONIAN_ERROR_HA:.3f} Ha cap
    ({error_2 / MAX_HAMILTONIAN_ERROR_HA * 100:.0f}%), and <code>E_gap</code> at {e_gap:.4f} Ha
    against a strict {MIN_E_GAP_HA} Ha floor (+{(e_gap / MIN_E_GAP_HA - 1) * 100:.1f}%). How far a
    candidate can move before one of them fires has not been measured.</p>
  </section>

  <section>
    <h2>Cost anatomy</h2>
    <p class="finding">A qubit saved is worth twice as much on <code>dqi</code>, whose cost carries
    Q²; the other two are linear in Q. That exponent, not the raw counts, decides which lever
    pays — and the three tasks span {toffoli_span:,.0f}× in Toffoli count, so they are not
    interchangeable targets.</p>
    <div class="tw"><table>
      <thead><tr><th>Task</th><th>Reference cost</th><th>Scored as</th><th>Q</th>
      <th>Toffolis/shot</th><th>Leading factor</th><th>Shots</th></tr></thead>
      <tbody>{reference_rows(reports)}</tbody>
    </table></div>
    <p class="note">Each row&rsquo;s cost is reconstructable from the columns beside it. For
    <code>invsqrt</code> the leading factor C = {float(invsqrt.detail["toffoli_equivalent"]):,.2f}
    already contains the Toffolis, so its cost is Q × C and the Toffoli column is shown for
    comparison only.</p>
  </section>

  <section>
    <h2>What the local loop costs</h2>
    <p class="finding">A full three-task verification takes about
    {sum(HAND_RECORDED["minutes"].values()):.0f} minutes on 10 CPU cores and is free and
    unlimited — against {HAND_RECORDED["submissions_per_day"]} rationed submissions a day. That
    ratio is the whole reason the dev set is worth building.</p>
    <div class="funnel">{run_log_rows(reports)}</div>
    <p class="note">Not a funnel yet: nothing has been eliminated, because no candidates have been
    generated. The per-stage denominators this section will carry only exist once a generator is
    running.</p>
  </section>

  <section>
    <h2>Checking the blend arithmetic</h2>
    <p class="finding">Feeding the public frontier&rsquo;s self-reported per-task costs through our
    scoring path reproduces the blend they reported, {PUBLIC_BLEND}, to
    {blend_delta_ppm:.1f} parts per million. That checks the ratio-and-geometric-mean arithmetic
    against a third party&rsquo;s independent calculation — it does <em>not</em> check the three
    per-task cost formulas, since both sides use the same self-reported costs.</p>
    <div class="tw"><table>
      <thead><tr><th>Task</th><th>Public best cost (third-party)</th>
      <th>Ratio vs our reference</th></tr></thead>
      <tbody>{ratio_rows(reports, calibration)}</tbody>
    </table></div>
    <p class="note">Computed with <code>cost.task_ratio</code> and <code>cost.score</code> against
    the measured calibration, giving {blend:.5f}. The agreement cannot be read as finer than the
    inputs: the rounding in those
    costs as quoted moves the blend by about ±{blend_floor_ppm:.1f} ppm on its own, so
    {blend_delta_ppm:.1f} ppm means &ldquo;agrees to the precision available&rdquo;, not a
    measured resolution.</p>
  </section>

  <section>
    <h2>Verdict ledger</h2>
    {verdict_blocks(verdicts)}
    <h3>Promotion gates</h3>
    <div class="tw"><table><tbody>{gate_rows(verdicts)}</tbody></table></div>
  </section>

  <section>
    <h2>Not established</h2>
    <ul>
      <li><strong>No compression has been attempted.</strong> Everything here is measurement
      apparatus.</li>
      <li><strong>Run-to-run noise is unmeasured.</strong> Every verifier seeds itself from the
      candidate&rsquo;s own bytes — <code>int(circuit_sha256, 16)</code>, or
      <code>sha256(hash + salt)</code> with a salt — so a re-run of one circuit is bit-identical
      and two candidates never share a draw. Comparing candidates needs an explicit
      <code>--seed</code>; sizing the across-seed noise is gate G1 and the next thing to run.</li>
      <li><strong>Every third-party number here is unverified by us.</strong> The per-task
      frontier costs and the {PUBLIC_BLEND} blend are self-reported in Kaggle topic 754174, and
      the blend check above rests entirely on them — if they are wrong, that check says nothing.
      The leaderboard scores and the {teams:,} team count are the competition&rsquo;s own.</li>
      <li><strong>How much <code>dfthcsos</code> headroom exists is unmeasured.</strong> The
      reference sits at 93% of the ANGLE qubit cap and 93% of the <code>error_2</code> budget, and
      nobody has measured how far a candidate can move before one rejects.</li>
      <li><strong>One leaderboard reading, on an input designed to score 1.0.</strong> Agreement
      there does not prove agreement on a compressed candidate.</li>
      <li><strong>Which task to commit to is open.</strong> The structural facts above bear on it;
      the decision does not follow from them.</li>
      <li><strong>Why {escaped} teams cleared 25 while half the field sits near
      {median:.1f}</strong> is unexplained, and which task holds that margin is the most valuable
      unknown in the competition.</li>
    </ul>
  </section>

  <footer>
    Built by <code>scripts/build_artifact.py</code> from the three official verifier outputs
    committed under <code>tests/data/</code>, the verdict ledger and leaderboard snapshot in
    <code>docs/experiments/</code>, and the formulas and constants imported from
    <code>quantum_compilation.cost</code>. The three wall-clock run times and the
    {HAND_RECORDED["submissions_per_day"]}-a-day submission budget are recorded by hand and live
    in <code>HAND_RECORDED</code>; third-party figures are labelled where they appear. Status
    {STATUS_DATE}.
  </footer>
</main>
</html>
"""


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(build())
    print(f"wrote {OUTPUT_PATH} ({OUTPUT_PATH.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
