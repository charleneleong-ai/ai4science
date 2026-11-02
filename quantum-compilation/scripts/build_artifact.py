"""Render the VQCC progress views from the verifiers' own outputs.

Numbers come from four committed files -- the official verifier runs over
`sample_submission/` in `tests/data/reference_<task>.json`, the verdict ledger in
`docs/experiments/verdicts.yaml`, the dated public-leaderboard snapshot in
`docs/experiments/leaderboard-<date>.json`, and the submitted parts' provenance record --
plus the formulas and `DEFAULT_CALIBRATION` imported from `quantum_compilation.cost`, so the
page cannot drift from the module the project scores with. `HAND_RECORDED` holds the handful
of figures that exist in no file (wall-clock run times, the competition's submission budget);
the footer says so.

    uv run --extra tracking python scripts/build_artifact.py

Writes `docs/artifact/index.html`, committed on every rebuild and republished as a Claude
Artifact at each significant result update, and `experiments/progress/vqcc/progress.png`,
the PR's visual aid, rendered by `autoresearch.render` from the same ledger.
"""

from __future__ import annotations

import html
import json
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, ClassVar

import yaml
from autoresearch.render import render

from quantum_compilation.cost import (
    AngleFloorSweep,
    clopper_pearson_upper_bound,
    BASELINE_FACTORIZATION,
    DEFAULT_CALIBRATION,
    ESCAPE_SCORE,
    FEMOCO_C0,
    MAX_ANGLE_QUBITS,
    MAX_DQI_FAILURE_RATE,
    NUM_ANGLE_WORDS,
    TASKS,
    ShapeResources,
    blend_factor,
    cost_cut_for_blend,
    dfthcsos_constraint_margins,
    dqi_repetition_factor,
    marginal_blend_per_unit,
    lambda_eff_floor,
    lcu_vectors,
    min_lcu_lambda,
    rotation_t_count,
    score,
    task_ratio,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = PROJECT_ROOT / "tests" / "data"
EXPERIMENTS_DIR = PROJECT_ROOT / "docs" / "experiments"
VERDICTS_PATH = EXPERIMENTS_DIR / "verdicts.yaml"
LEADERBOARD_PATH = EXPERIMENTS_DIR / "leaderboard-2026-10-31.json"
PARTS_PATH = EXPERIMENTS_DIR / "57834335-parts.yaml"
SEED_RECORD_PATH = EXPERIMENTS_DIR / "57834335-seeds.yaml"
OUTPUT_PATH = PROJECT_ROOT / "docs" / "artifact" / "index.html"
PROGRESS_PATH = PROJECT_ROOT / "experiments" / "progress" / "vqcc" / "progress.png"

STATUS_DATE = "2026-11-02"
# A bar within this fraction of its limit counts as pinned, for the chart's colouring and for
# the count the prose quotes. One constant so the two cannot disagree.
PINNED_FRACTION = 0.999


# Wall-clock minutes per verifier on 10 CPU cores, and the competition's daily submission
# budget. Neither is in any result file; the footer declares them as recorded by hand.
HAND_RECORDED = {"minutes": {"invsqrt": 5.7, "dqi": 2.4, "dfthcsos": 1.0}, "submissions_per_day": 5}

# Read off the submitted ANGLE circuits and off two rebuilds of the same `angle_words`, none of
# which any result file records: the gate histogram of `angle_k_13.kmb` and the register each
# of its CX gates targets, and the Toffolis per pass that kickmix's own `init_lookup` and the
# starter kit's segment-tree QROM cost for that same table. Declared here rather than written
# into the prose, so the page, the README and the `r3-crack-the-angle-loaders` verdict row
# quote one source; every ratio between them is computed.
LOADER_REBUILDS = {
    "init_lookup_toffoli": 30_690,
    "qrom_toffoli": 18_270,
    "circuit_cx": 33_540,
    "selection_cx": 27_444,
    "cheapest_circuit_toffoli": 258,
    "dearest_circuit_toffoli": 301,
    "extra_workspace_qubits": 65,
}

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

# How each task's cost is assembled, read against cost.py's three cost functions: the written
# form, the reported field carrying its leading factor, and that factor's symbol.
SCORED_AS = {
    "invsqrt": ("Q × C", "toffoli_equivalent", "C"),
    "dqi": ("R(ε₊) × Q² × T", "repetition_factor", "R"),
    "dfthcsos": ("λ_eff × T × Q", "lambda_eff", "λ_eff"),
}

QUANTITY_LABELS = {
    "hamiltonian_error_ha": "Hamiltonian error",
    "e_gap_ha": "E_gap",
    "angle_qubits": "widest ANGLE circuit",
    "lcu_lambda": "one-norm Λ",
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


@dataclass(frozen=True)
class Standing:
    """Where one score sits in a leaderboard snapshot, its tie block included.

    Written once and used twice -- for our score and for the 1.00000 block the reference
    submission shares -- so the page reports the move with the same honesty at both ends.
    """

    # Leaderboard scores are published rounded, so an exact float comparison would split a
    # tie block; this is half of the last published digit.
    TOLERANCE: ClassVar[float] = 5e-6

    scores: list[float]
    score: float

    @property
    def teams(self) -> int:
        return len(self.scores)

    @property
    def above(self) -> int:
        return sum(1 for s in self.scores if s > self.score + self.TOLERANCE)

    @property
    def tied(self) -> int:
        return sum(1 for s in self.scores if abs(s - self.score) < self.TOLERANCE)

    @property
    def best(self) -> int:
        return self.above + 1

    @property
    def worst(self) -> int:
        return self.above + self.tied

    @property
    def leader(self) -> float:
        return self.scores[0]


def scored_submissions(verdicts: dict[str, Any]) -> list[dict[str, Any]]:
    """The ledger entries carrying a leaderboard score, oldest first.

    `verdicts.yaml` is the single source for the page's headline, the progress plot and the
    W&B report, so each score is recorded once beside the experiment that earned it.
    """
    return [e for e in verdicts.get("experiments", []) if e.get("score") is not None]


def write_progress_png(verdicts: dict[str, Any]) -> Path:
    """The PR's visual aid: each submission's leaderboard score, with the running best.

    `autoresearch.render` owns the chart, so all this does is map the ledger onto the rows it
    reads. `status` is what it draws the running-best step line from, and every submission so
    far was promoted, so the first is the BASELINE and the rest are KEEP.
    """
    rows = [
        {
            "experiment": index,
            "score": float(entry["score"]),
            "status": "BASELINE" if index == 0 else "KEEP",
            "name": entry["name"],
        }
        for index, entry in enumerate(scored_submissions(verdicts))
    ]
    return render(
        rows=rows,
        out=PROGRESS_PATH,
        title="VQCC leaderboard score per submission",
        score_label="Leaderboard score (higher is better)",
        score_fmt="{:.5f}",
        label_field="name",
        y_margin=0.12,
    )


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


def is_pinned(margin: dict[str, Any]) -> bool:
    return margin["fraction_of_limit_used"] > PINNED_FRACTION


def binding_margin(key: str, margin: dict[str, Any]) -> dict[str, Any]:
    """`margin`, with `lcu_lambda` re-based onto the floor that actually binds it.

    The constraint table carries `MIN_LAMBDA = 49.0` because that is the verifier's rejection
    rule, and `dfthcsos_rejections` must mirror the verifier. But it never binds: `cost`'s
    `min_lcu_lambda` is 281.52 for this reference, so plotting against 49.0 draws the one-norm
    with 84% headroom when it has 6.1% -- the least of the four, not the most.
    """
    if key != "lcu_lambda":
        return margin
    floor = min_lcu_lambda(FEMOCO_C0)
    return margin | {
        "limit": floor,
        "fraction_of_limit_used": floor / margin["value"],
        "note": f"derived; the {margin['limit']:g} rule never binds",
    }


def constraint_label(key: str, margin: dict[str, Any]) -> str:
    """"E_gap vs 3.914898 floor" -- the limit and its direction read off the record itself."""
    word = "cap" if margin["sense"] == "at_or_below" else "floor"
    limit = margin["limit"]
    note = f" ({margin['note']})" if "note" in margin else ""
    return f"{QUANTITY_LABELS[key]} vs {limit:.6g} {word}{note}"


def svg(width: int, height: int, aria: str, body: str) -> str:
    return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{aria}">{body}</svg>'


def hgrid(x0: float, x1: float, y: float, label: str) -> str:
    return (
        f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" class="grid"/>'
        f'<text x="{x0 - 10:.1f}" y="{y + 6:.1f}" text-anchor="end" class="ct">{label}</text>'
    )


def field_svg(standing: Standing) -> str:
    """The whole field, sorted: the plateau, the cliff, and which block we are in."""
    scores = standing.scores
    width, height = 680, 400
    pad_l, pad_r, pad_t, pad_b = 74, 26, 24, 62
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b
    teams = standing.teams
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
    x_best, x_worst = place(standing.best, 0)[0], place(standing.worst, 0)[0]
    our_y = place(standing.best, standing.score)[1]
    leader_x, leader_y = place(1, standing.leader)
    return svg(
        width,
        height,
        f"All {teams:,} public leaderboard scores in rank order: a dense plateau near 23, a "
        f"cliff, and our score of {standing.score:.2f} on the plateau's top edge",
        grid
        + f'<rect x="{x_best:.1f}" y="{pad_t}" width="{max(3.0, x_worst - x_best):.1f}" '
        f'height="{plot_h}" fill="var(--warn)" fill-opacity="0.18"/>'
        + ticks
        + f'<path d="{" ".join(path)}" fill="none" stroke="var(--accent)" stroke-width="3"/>'
        + f'<circle cx="{leader_x:.1f}" cy="{leader_y:.1f}" r="7" fill="var(--accent)"/>'
        + f'<text x="{leader_x + 14:.1f}" y="{leader_y + 7:.1f}" class="cv">'
        f"leader {standing.leader:.2f}</text>"
        + f'<text x="{x_worst + 12:.1f}" y="{our_y - 16:.1f}" class="cv">'
        f"us {standing.score:.2f}</text>"
        + f'<text x="{x_worst + 12:.1f}" y="{our_y + 4:.1f}" class="ct">'
        f"rank {standing.best:,} of {teams:,}</text>"
        + f'<text x="{pad_l + plot_w / 2:.0f}" y="{height - 8}" text-anchor="middle" '
        f'class="cax">rank, of {teams:,} teams</text>'
        + f'<text x="22" y="{pad_t + plot_h / 2:.0f}" class="cax" '
        f'transform="rotate(-90 22 {pad_t + plot_h / 2:.0f})" text-anchor="middle">score</text>',
    )


def headroom_svg(margins: dict[str, dict[str, Any]]) -> str:
    """How much of each hard constraint the plateau's chemistry parameterisation already spends.

    Two bars sit on their limit, a third is two qubits from its cap, and the fourth -- the
    one-norm, which looked like the exception -- is at 93.9% of the floor that actually binds
    it. Nothing is left to tune inside this factorization, and `r1-price-the-factorization-
    refit` showed a different one cannot pay either: the gain has to come from `T` and `Q`.
    """
    margins = {key: binding_margin(key, margin) for key, margin in margins.items()}
    order = sorted(margins, key=lambda key: -margins[key]["fraction_of_limit_used"])
    width = 680
    pad_l, pad_r, pad_t, row_h = 300, 86, 30, 46
    height = pad_t + row_h * len(order) + 44
    plot_w = width - pad_l - pad_r

    body = [
        f'<line x1="{pad_l + plot_w:.1f}" y1="{pad_t - 8}" x2="{pad_l + plot_w:.1f}" '
        f'y2="{pad_t + row_h * len(order):.1f}" class="grid"/>',
        f'<text x="{pad_l + plot_w:.1f}" y="{height - 22}" text-anchor="middle" '
        f'class="ct">limit</text>',
        f'<text x="{pad_l:.1f}" y="{height - 22}" text-anchor="middle" class="ct">0</text>',
    ]
    for index, key in enumerate(order):
        used = margins[key]["fraction_of_limit_used"]
        y = pad_t + index * row_h
        colour = "var(--warn)" if is_pinned(margins[key]) else "var(--accent)"
        body += [
            f'<text x="{pad_l - 14:.1f}" y="{y + 20:.1f}" text-anchor="end" '
            f'class="ct">{constraint_label(key, margins[key])}</text>',
            f'<rect x="{pad_l:.1f}" y="{y + 6:.1f}" width="{plot_w:.1f}" height="22" '
            f'fill="var(--grid)" fill-opacity="0.35"/>',
            f'<rect x="{pad_l:.1f}" y="{y + 6:.1f}" width="{max(2.0, used * plot_w):.1f}" '
            f'height="22" fill="{colour}"/>',
            f'<text x="{pad_l + plot_w + 10:.1f}" y="{y + 22:.1f}" class="cv">'
            f'{used * 100:.3f}%</text>',
        ]
    pinned = sum(1 for key in order if is_pinned(margins[key]))
    return svg(
        width,
        height,
        f"{len(order)} hard constraints on the chemistry task, each as a fraction of its "
        f"limit; {pinned} are spent to more than "
        f"{(1 - PINNED_FRACTION) * 100:.1f}% of their limit",
        "".join(body),
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


def toffoli_split_svg(entries: Sequence[tuple[str, float, float]]) -> str:
    """One stacked bar per submission: the ANGLE Toffolis, then the heuristic ones.

    Both on a shared scale, because the comparison that matters is absolute -- the plateau did
    not merely shift the ratio, it cut the total by 5.7x. The ANGLE segment is the accented one
    because it is the only part a submission can lower without spending accuracy.
    """
    width, bar_h, gap, pad_l, pad_r = 680, 46, 52, 150, 150
    plot = width - pad_l - pad_r
    scale = plot / max(free + bought for _, free, bought in entries)
    rows = []
    for index, (label, free, bought) in enumerate(entries):
        y = index * (bar_h + gap) + 20
        w_free, w_bought = free * scale, bought * scale
        rows.append(
            f'<text x="{pad_l - 14}" y="{y + 28}" text-anchor="end" class="cl">{label}</text>'
            f'<rect x="{pad_l}" y="{y}" width="{w_free:.1f}" height="{bar_h}" rx="4" '
            f'fill="var(--accent)"/>'
            f'<rect x="{pad_l + w_free:.1f}" y="{y}" width="{w_bought:.1f}" height="{bar_h}" '
            f'rx="4" fill="var(--neutral-mark)" fill-opacity="0.45"/>'
            f'<text x="{pad_l + 12}" y="{y + 29}" class="cvi">{free / (free + bought):.0%}</text>'
            f'<text x="{pad_l + w_free + w_bought + 12:.1f}" y="{y + 21}" class="cv">'
            f"{free + bought:,.0f}</text>"
            f'<text x="{pad_l + w_free + w_bought + 12:.1f}" y="{y + 39}" class="ct">'
            f"{free:,.0f} free + {bought:,.0f} bought</text>"
        )
    legend = (
        f'<text x="{pad_l}" y="{len(entries) * (bar_h + gap) + 14}" class="ct">'
        "accented: ANGLE circuits, which cost no accuracy &#183; grey: fixed by the shape and brot"
        "</text>"
    )
    return svg(
        width,
        len(entries) * (bar_h + gap) + 24,
        "Half of the plateau submission's Toffoli count is ANGLE circuits, which cost no "
        "accuracy budget",
        "".join(rows) + legend,
    )


def angle_floor_svg(sweep: AngleFloorSweep) -> str:
    """Blend against batch size, at the AND counting floor and at the submitted AND rate.

    Two lines on one scale because the gap between them is the finding: the floor clears the
    escape band everywhere the contract allows and the measured construction clears none of
    it, so what separates them is circuit quality rather than budget. The window is taken from
    the data, including the two reference rules, so nothing plotted can fall outside it.
    """
    width, height = 680, 360
    pad_l, pad_r, pad_t, pad_b = 66, 26, 26, 62
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b
    plotted = [v for row in sweep.rows for v in (row.at_floor, row.at_rate)]
    y_lo = math.floor(min(*plotted, sweep.blend, sweep.target) * 2) / 2
    y_hi = math.ceil(max(*plotted, sweep.blend, sweep.target) * 2) / 2
    x_lo, x_hi = sweep.rows[0].batch_size, sweep.rows[-1].batch_size

    def place(batch: int, value: float) -> tuple[float, float]:
        return (
            pad_l + (batch - x_lo) / (x_hi - x_lo) * plot_w,
            pad_t + plot_h - (value - y_lo) / (y_hi - y_lo) * plot_h,
        )

    def line(points: Sequence[tuple[float, float]], colour: str) -> str:
        path = " ".join(
            f"{'M' if step == 0 else 'L'}{x:.1f} {y:.1f}" for step, (x, y) in enumerate(points)
        )
        dots = "".join(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colour}"/>' for x, y in points
        )
        return f'<path d="{path}" fill="none" stroke="{colour}" stroke-width="3"/>{dots}'

    x_right = pad_l + plot_w
    target_y, blend_y = place(x_lo, sweep.target)[1], place(x_lo, sweep.blend)[1]
    best = sweep.best_floor
    bx, by = place(best.batch_size, best.at_floor)
    grid = "".join(
        hgrid(pad_l, x_right, place(x_lo, tick)[1], f"{tick:.0f}")
        for tick in range(math.ceil(y_lo), math.floor(y_hi) + 1)
    )
    ticks = "".join(
        f'<text x="{place(row.batch_size, y_lo)[0]:.1f}" y="{height - 34}" '
        f'text-anchor="middle" class="ct">{row.batch_size}</text>'
        for row in sweep.rows
    )
    return svg(
        width,
        height,
        "Blend against the ANGLE batch size, at the AND counting floor and at the submitted "
        "rate, with the escape score marked",
        grid
        + ticks
        + f'<line x1="{pad_l}" y1="{target_y:.1f}" x2="{x_right}" y2="{target_y:.1f}" '
        'class="pole"/>'
        + f'<text x="{x_right - 6:.1f}" y="{target_y + 19:.1f}" text-anchor="end" class="cl2">'
        f"escape band, {sweep.target:g}</text>"
        + f'<line x1="{pad_l}" y1="{blend_y:.1f}" x2="{x_right}" y2="{blend_y:.1f}" '
        'class="marker"/>'
        + f'<text x="{x_right - 6:.1f}" y="{blend_y + 18:.1f}" text-anchor="end" class="ct">'
        f"submitted, {sweep.blend:.2f}</text>"
        + line([place(row.batch_size, row.at_floor) for row in sweep.rows], "var(--accent)")
        + line([place(row.batch_size, row.at_rate) for row in sweep.rows], "var(--warn)")
        + f'<text x="{bx + 10:.1f}" y="{by - 8:.1f}" class="cv">'
        f"J={best.batch_size}, {best.at_floor:.2f}</text>"
        + f'<text x="{x_right - 6:.1f}" y="{pad_t + plot_h - 9:.1f}" text-anchor="end" '
        'class="ct">accent: counting floor &#183; warm: submitted rate</text>'
        + f'<text x="{pad_l + plot_w / 2:.0f}" y="{height - 8}" text-anchor="middle" '
        'class="cax">J, angle-word registers per circuit</text>',
    )


def angle_share_svg(report: TaskReport) -> str:
    """Two true partitions: what the ANGLE circuits are of the Toffolis, and of the width.

    Only additive decompositions go in a part-of-whole bar. The scored total width (730) is
    *not* one -- `temporary_qubits` is a max over the three families, not a sum -- so the
    bars show the Toffoli partition and the persistent-qubit partition, and the prose below
    carries the max.
    """
    res = ShapeResources.from_report(report.detail)
    bars = [
        ("Toffolis", res.angle_toffoli, res.toffoli, "of all Toffolis"),
        (
            "Persistent qubits",
            float(res.angle_word_qubits),
            float(res.heuristic_persistent_qubits + res.angle_word_qubits),
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


# Not an f-string, so the braces are the CSS's own and need no doubling.
STYLE = """\
/* Instrument readout: one measured column, data in mono, the finding above each chart. */
:root {
  --bg: #f7f8f7; --surface: #ffffff; --fg: #15201c; --muted: #5b6b64;
  --line: #dfe4e1; --grid: #b4c0ba; --accent: #0f6f5c; --accent-soft: #e4f0ec;
  --warn: #b5700f; --warn-soft: #fbeedd; --bad: #a3342a; --bad-soft: #f8e5e2;
  --neutral-mark: #8a9890;
  --display: "IBM Plex Serif", Georgia, serif;
  --body: "IBM Plex Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0e1412; --surface: #151d1a; --fg: #e6ece9; --muted: #93a39c;
    --line: #26322e; --grid: #46554f; --accent: #4fc4a6; --accent-soft: #16302a;
    --warn: #e0a449; --warn-soft: #33270f; --bad: #e07a6e; --bad-soft: #351e1b;
    --neutral-mark: #5d6d67;
    color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #0e1412; --surface: #151d1a; --fg: #e6ece9; --muted: #93a39c;
  --line: #26322e; --grid: #46554f; --accent: #4fc4a6; --accent-soft: #16302a;
  --warn: #e0a449; --warn-soft: #33270f; --bad: #e07a6e; --bad-soft: #351e1b;
  --neutral-mark: #5d6d67;
  color-scheme: dark;
}
* { box-sizing: border-box; }
body {
  background: var(--bg); color: var(--fg); font-family: var(--body);
  font-size: 15px; line-height: 1.6; margin: 0;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 40px 72px;
}
main { max-width: 760px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }
h1 { font-family: var(--display); font-size: clamp(27px, 5vw, 38px); line-height: 1.15;
      margin: 0 0 12px; text-wrap: balance; font-weight: 600; }
h2 { font-family: var(--display); font-size: 21px; margin: 0 0 4px; font-weight: 600; }
h3 { font-size: 13px; text-transform: uppercase; letter-spacing: 0.08em;
      color: var(--muted); margin: 0 0 10px; font-weight: 600; }
p { margin: 0 0 12px; }
section { display: flex; flex-direction: column; gap: 14px; }
.status { font-family: var(--mono); font-size: 12.5px; color: var(--muted); }
.lede { font-size: 18px; line-height: 1.5; }
.strip { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; }
.tile { background: var(--surface); border: 1px solid var(--line); border-radius: 8px;
         padding: 14px 16px; display: flex; flex-direction: column; gap: 3px; min-width: 0; }
.tile .k { font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.07em;
            color: var(--muted); font-weight: 600; }
.tile .v { font-family: var(--mono); font-size: 23px; font-weight: 600;
            font-variant-numeric: tabular-nums; }
.tile .n { font-size: 12.5px; color: var(--muted); }
.tile.flag { border-color: var(--warn); background: var(--warn-soft); }
.finding { font-size: 17px; line-height: 1.45; font-weight: 500; margin: 0;
            border-left: 3px solid var(--accent); padding: 2px 0 2px 14px; }
figure { margin: 0; background: var(--surface); border: 1px solid var(--line);
          border-radius: 8px; padding: 16px; display: flex; flex-direction: column; gap: 10px; }
figure svg { width: 100%; height: auto; max-width: 100%; display: block; }
figcaption { font-size: 12.5px; color: var(--muted); }
.grid { stroke: var(--grid); stroke-width: 1; }
.pole { stroke: var(--bad); stroke-width: 2; stroke-dasharray: 5 4; }
.marker { stroke: var(--warn); stroke-width: 1.5; stroke-dasharray: 4 4; }
text { font-family: var(--body); fill: var(--fg); }
.cl { font-size: 20px; fill: var(--fg); }
.cl2 { font-size: 19px; fill: var(--muted); font-weight: 600; }
.ct { font-size: 18px; fill: var(--muted); font-family: var(--mono); }
.cax { font-size: 19px; fill: var(--muted); }
.cv { font-size: 20px; fill: var(--fg); font-family: var(--mono); font-weight: 600; }
.cvi { font-size: 20px; fill: var(--bg); font-family: var(--mono); font-weight: 600; }
.tw { overflow-x: auto; border: 1px solid var(--line); border-radius: 8px;
       background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { text-align: left; padding: 9px 12px; border-bottom: 1px solid var(--line);
          vertical-align: top; }
th { font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.06em;
      color: var(--muted); font-weight: 600; white-space: nowrap; }
tr:last-child td { border-bottom: none; }
td.num { font-family: var(--mono); font-variant-numeric: tabular-nums; white-space: nowrap; }
code { font-family: var(--mono); font-size: 0.92em; }
.pill { display: inline-block; font-family: var(--mono); font-size: 11px; font-weight: 600;
         padding: 2px 7px; border-radius: 4px; letter-spacing: 0.04em; }
.pill.ok { background: var(--accent-soft); color: var(--accent); }
.pill.warn { background: var(--warn-soft); color: var(--warn); }
.pill.bad { background: var(--bad-soft); color: var(--bad); }
.pill.run { background: var(--surface); color: var(--muted); border: 1px solid var(--line); }
.funnel { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; }
.fcell { border: 1px solid var(--line); border-left: 3px solid var(--neutral-mark);
          border-radius: 6px; padding: 12px 14px; background: var(--surface);
          display: flex; flex-direction: column; gap: 2px; min-width: 0; }
.ftask { font-family: var(--mono); font-weight: 600; font-size: 14px; }
.fpass { font-family: var(--mono); font-size: 11px; font-weight: 600; color: var(--muted);
          letter-spacing: 0.05em; }
.fmeta { font-size: 12px; color: var(--muted); font-family: var(--mono); }
.verdict { border: 1px solid var(--line); border-radius: 8px; background: var(--surface);
            padding: 14px 16px; }
.verdict header { display: flex; gap: 10px; align-items: baseline; flex-wrap: wrap; }
.verdict header code { font-size: 13.5px; font-weight: 600; }
.verdict p { margin: 0; font-size: 14px; }
.verdict .vmeta { font-size: 12px; color: var(--muted); font-family: var(--mono);
                   margin: 2px 0 8px; }
ul { margin: 0; padding-left: 20px; }
li { margin-bottom: 8px; }
li::marker { color: var(--neutral-mark); }
.note { border-left: 3px solid var(--neutral-mark); padding: 2px 0 2px 14px;
         color: var(--muted); font-size: 14px; }
footer { border-top: 1px solid var(--line); padding-top: 16px; font-size: 12.5px;
          color: var(--muted); }
a { color: var(--accent); }
"""


@dataclass(frozen=True)
class Artifact:
    """The page, and every figure on it, derived from the committed result files.

    One object so each section reads the figures it needs off `self` instead of thirty
    locals being threaded through a single template, and so a figure used by two sections --
    the rank, the local blend, the chemistry margins -- is computed once in one place.
    """

    reports: dict[str, TaskReport]
    verdicts: dict[str, Any]
    board: dict[str, Any]
    parts: dict[str, Any]
    seed_runs: dict[str, Any]

    @classmethod
    def load(cls) -> Artifact:
        board = json.loads(LEADERBOARD_PATH.read_text())
        # The page calls the team count the competition's own, and reads it as the number of
        # scores in the snapshot, so the snapshot's own figure has to agree with that.
        if board["teams"] != len(board["scores"]):
            raise ValueError(
                f"{LEADERBOARD_PATH.name} reports {board['teams']} teams but carries "
                f"{len(board['scores'])} scores"
            )
        return cls(
            reports=load_reports(),
            verdicts=yaml.safe_load(VERDICTS_PATH.read_text()),
            board=board,
            parts=yaml.safe_load(PARTS_PATH.read_text()),
            seed_runs=yaml.safe_load(SEED_RECORD_PATH.read_text()),
        )

    # Only figures two or more sections share are properties; anything one section uses is a
    # local in that section, where it reads next to the prose quoting it. Every score on the
    # page goes through `cost.DEFAULT_CALIBRATION`, the module's own measurement of the same
    # three verifier outputs loaded here, so the page cannot score against a calibration the
    # project does not use -- `tests/test_verify.py` pins the two to each other.

    @property
    def scores(self) -> list[float]:
        return self.board["scores"]

    @property
    def our_score(self) -> float:
        """The latest leaderboard reading -- the one figure on the page not derived locally.

        Recorded once in `verdicts.yaml` beside the submission that earned it, rather than as
        a constant here; the local blend it is compared against is recomputed, so the page's
        central claim checks itself instead of asserting a transcription.
        """
        return float(scored_submissions(self.verdicts)[-1]["score"])

    @cached_property
    def standing(self) -> Standing:
        return Standing(self.scores, self.our_score)

    @property
    def median(self) -> float:
        return self.scores[self.standing.teams // 2]

    @property
    def plateau_hi(self) -> float:
        """The 10th-percentile score: the top of the band half the field is packed into."""
        return self.scores[int(self.standing.teams * 0.10)]

    @property
    def escaped(self) -> int:
        return sum(1 for s in self.scores if s >= ESCAPE_SCORE)

    @property
    def part_costs(self) -> dict[str, float]:
        """What our verifier measured for the circuits we actually submitted.

        Distinct from `reports`, whose costs come from `sample_submission/`.
        """
        return {task: self.parts["parts"][task]["our_measurement"]["cost"] for task in TASKS}

    @cached_property
    def seed_spread(self) -> dict[str, dict[str, float]]:
        """How far each task's cost moved across the seeds, recomputed from the record.

        The record stores only what the verifiers reported, so these are derived here rather
        than read: a stored spread could disagree with the rows it came from.
        """
        seeds = self.seed_runs["seeds"]
        series = {task: [s["costs"][task] for s in seeds] for task in TASKS}
        series["blend"] = [score(s["costs"], DEFAULT_CALIBRATION) for s in seeds]
        series["toffolis"] = [s["dqi_detail"]["toffoli_count"] for s in seeds]
        return {
            label: {
                "min": min(values),
                "max": max(values),
                "range": max(values) - min(values),
                "ppm": (max(values) - min(values)) / (sum(values) / len(values)) * 1e6,
            }
            for label, values in series.items()
        }

    @cached_property
    def failed_shots(self) -> list[int]:
        """Shots the dqi decoder dropped, per seed -- the whole of that task's variance."""
        shots = self.seed_runs["shots_per_run"]
        return [round(s["dqi_detail"]["empirical_error_rate"] * shots) for s in self.seed_runs["seeds"]]

    @cached_property
    def failure_interval(self) -> tuple[float, float]:
        """95% Clopper-Pearson interval on how often a run drops at least one shot.

        The same distribution the dqi verifier scores on, so it comes from the same helper:
        a two-sided interval is its one-sided bound read from both ends.
        """
        runs = len(self.failed_shots)
        hits = sum(1 for f in self.failed_shots if f)
        tail = math.log2(0.025)
        return (
            1.0 - clopper_pearson_upper_bound(runs - hits, runs, log2_delta=tail),
            clopper_pearson_upper_bound(hits, runs, log2_delta=tail),
        )

    @cached_property
    def local_blend(self) -> float:
        """The whole headline claim is that this predicted the leaderboard, so it is
        recomputed from the submitted parts' measured costs rather than transcribed."""
        return score(self.part_costs, DEFAULT_CALIBRATION)

    # -- rows and blocks the tables are made of --

    def reference_rows(self) -> str:
        rows = []
        for task in TASKS:
            report = self.reports[task]
            formula, field_name, symbol = SCORED_AS[task]
            value = float(report.detail[field_name])
            rows.append(
                f"<tr><td><code>{task}</code></td>"
                f'<td class="num">{fmt(report.cost)}</td>'
                f"<td><code>{formula}</code></td>"
                f'<td class="num">{int(report.detail["num_qubits"]):,}</td>'
                f'<td class="num">{report.toffolis:,.1f}</td>'
                f'<td class="num">{symbol} = {value:,.4f}</td>'
                f'<td class="num">{report.shots:,}</td></tr>'
            )
        return "".join(rows)

    def ratio_rows(self) -> str:
        return "".join(
            f"<tr><td><code>{task}</code></td>"
            f'<td class="num">{PUBLIC_COST_QUOTED[task]}</td>'
            f'<td class="num">{task_ratio(task, PUBLIC_COSTS[task], DEFAULT_CALIBRATION):,.3f}×'
            f"</td></tr>"
            for task in TASKS
        )

    def seed_rows(self) -> str:
        """Costs render through `fmt`; the blend is a score, so it keeps the decimals its
        whole spread lives in."""
        spread = self.seed_spread
        rows = []
        for label, reading in (
            ("dfthcsos", "identical every seed &mdash; deterministic"),
            ("invsqrt", "below 1 ppm &mdash; effectively deterministic"),
            ("dqi", "one discrete step: 0 or 1 failed shots"),
            ("blend", "carried entirely by <code>dqi</code>"),
        ):
            ppm = spread[label]["ppm"]
            lo, hi = spread[label]["min"], spread[label]["max"]
            shown = (f"{lo:,.5f}", f"{hi:,.5f}") if label == "blend" else (fmt(lo), fmt(hi))
            cell = (
                f'<span class="pill warn">{ppm:,.0f} ppm</span>'
                if ppm > 100
                else f"{ppm:,.2f} ppm"
            )
            rows.append(
                f"<tr><td><code>{label}</code></td>"
                f'<td class="num">{shown[0]}</td><td class="num">{shown[1]}</td>'
                f'<td class="num">{cell}</td><td>{reading}</td></tr>'
            )
        return "".join(rows)

    def run_log_rows(self) -> str:
        return "".join(
            f'<div class="fcell"><div class="ftask">{task}</div>'
            f'<div class="fpass">PASS</div>'
            f'<div class="fmeta">{self.reports[task].shots:,} shots</div>'
            f'<div class="fmeta">{HAND_RECORDED["minutes"][task]} min · '
            f"cost {fmt(self.reports[task].cost)}</div></div>"
            for task in TASKS
        )

    def verdict_blocks(self) -> str:
        tone = {"PROMOTE": "ok", "KEEP": "warn", "DISCARD": "bad", "RUNNING": "run"}
        blocks = []
        for item in self.verdicts.get("experiments", []):
            verdict = str(item.get("verdict", "RUNNING"))
            pr = item.get("pr")
            link = (
                f' · <a href="https://github.com/charleneleong-ai/ai4science/pull/{pr}">'
                f"PR #{pr}</a>"
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

    def gate_rows(self) -> str:
        return "".join(
            f"<tr><td><code>{html.escape(name)}</code></td><td>{inline_markup(text)}</td></tr>"
            for name, text in self.verdicts.get("gates", {}).items()
        )

    # -- the sections, in the order they are read --

    def header(self) -> str:
        return f"""
  <header>
    <p class="status">GOOGLE QUANTUM AI · VERIFIABLE QUANTUM CIRCUIT COMPILATION · {STATUS_DATE}</p>
    <h1>On the plateau</h1>
    <p class="lede">The cost model is read off the three shipped verifiers rather than the
    Evaluation page, and every formula reproduces the verifier&rsquo;s own score exactly. That
    local loop has now predicted the leaderboard twice — on the reference
    (1.00000&nbsp;→&nbsp;1.00000) and on an input it had never seen
    ({self.local_blend:.5f}&nbsp;→&nbsp;<strong>{self.our_score:.5f}</strong>, agreeing to every digit
    either is reported to) — so candidates can be ranked for free against a rationed
    {HAND_RECORDED["submissions_per_day"]} submissions a day. The score itself is assembled from
    public Apache-2.0 parts, not from any compression of our own; what it bought is a
    measurement of why the plateau is where it is.</p>
  </header>
"""

    def strip(self) -> str:
        standing = self.standing
        # The 1.00000 block, which 81 teams shared by submitting the reference unchanged.
        reference = Standing(self.scores, 1.0)
        gap_score = standing.leader / self.our_score
        return f"""
    <div class="strip">
      <div class="tile"><span class="k">Our score</span>
        <span class="v">{self.our_score:.5f}</span>
        <span class="n">local said {self.local_blend:.5f} · submission
        {self.parts["submission_ref"]}</span></div>
      <div class="tile"><span class="k">Rank</span>
        <span class="v">{standing.best:,}–{standing.worst:,}</span>
        <span class="n">of {standing.teams:,} · up from the {reference.tied}-team block at
        1.00000 ({reference.best:,}–{reference.worst:,})</span></div>
      <div class="tile flag"><span class="k">Gap to leader</span>
        <span class="v">×{gap_score:.2f}</span>
        <span class="n">≈{gap_score**3:,.1f}× on the product of the three cost ratios</span></div>
      <div class="tile flag"><span class="k">Compression attempted</span>
        <span class="v">none</span><span class="n">every circuit here is someone else&rsquo;s</span></div>
    </div>
    <p class="note">This is an assembly, not a result: the public parts market carries a better
    circuit than the reference for all three tasks, so the best per-task pick was fetched,
    SHA-256 checked against its manifest and scored locally before a submission was spent. It
    establishes the floor to beat and credits nine authors in
    <code>docs/experiments/{PARTS_PATH.name}</code>. The ×{gap_score:.2f} still owed to the
    leader is the whole remaining problem.</p>
"""

    def where_we_stand(self) -> str:
        standing = self.standing
        return f"""
    <h2>Where we stand</h2>
    <p class="finding">The field is bimodal. Half of it is packed between {self.median:.2f} and
    {self.plateau_hi:.2f} — a band {self.plateau_hi - self.median:.1f} wide — and only
    {self.escaped} teams have escaped past {ESCAPE_SCORE:g}.
    We now sit at rank {standing.best:,}–{standing.worst:,} of {standing.teams:,}, sharing our
    exact score with {standing.tied - 1} other teams, on the top edge of that band rather than
    below it.</p>
    <figure>
      {field_svg(standing)}
      <figcaption>Every public score on {STATUS_DATE}, in rank order, from the competition&rsquo;s
      own leaderboard snapshotted in
      <code>docs/experiments/leaderboard-{STATUS_DATE}.json</code>. The shaded band marks our
      position.</figcaption>
    </figure>
    <p>The crowded plateau is the useful part, and the section below says what it is made of: not
    a field that merely copied one notebook, but one that has run a single factorization up
    against every constraint the task has. The {self.escaped} teams above {ESCAPE_SCORE:g} found
    something that parameterisation cannot reach.</p>
    <p class="note">Why this assembly scores above the {PUBLIC_BLEND} blend quoted on
    {STATUS_DATE[:-3]}-30: the parts published since are cheaper on all three tasks —
    <code>invsqrt</code> {self.part_costs["invsqrt"]:,.0f} against
    {PUBLIC_COST_QUOTED["invsqrt"]}, <code>dqi</code> {fmt(self.part_costs["dqi"])} against
    {PUBLIC_COST_QUOTED["dqi"]}, <code>dfthcsos</code> {fmt(self.part_costs["dfthcsos"])} against
    {PUBLIC_COST_QUOTED["dfthcsos"]}. Their authors made them cheaper; we picked newer
    ones.</p>
"""

    def why_the_plateau_is_a_wall(self) -> str:
        plateau = self.parts["parts"]["dfthcsos"]
        margins = plateau["constraint_margins"]
        error, e_gap = margins["hamiltonian_error_ha"], margins["e_gap_ha"]
        angle, lcu = margins["angle_qubits"], margins["lcu_lambda"]
        e_gap_slack = (e_gap["value"] - e_gap["limit"]) * 1e6
        r, b, c = (plateau[key] for key in ("R", "B", "C"))
        br, bb, bc = BASELINE_FACTORIZATION
        vectors, ref_vectors = lcu_vectors(r, b), lcu_vectors(br, bb)
        capacity, ref_capacity = r * c, br * bc
        # The baseline's own margins, off the same table, for the accuracy it has not spent.
        reference = self.reports["dfthcsos"]
        ref_error = dfthcsos_constraint_margins(reference.detail)["hamiltonian_error_ha"][
            "fraction_of_limit_used"
        ]
        pinned = [QUANTITY_LABELS[key] for key in margins if is_pinned(margins[key])]
        return f"""
    <h2>Why the plateau is a wall</h2>
    <p class="finding">The chemistry circuits everyone is standing on have <strong>almost no
    headroom left against the constraints that reject outright</strong>: the accuracy term is
    at {error["fraction_of_limit_used"] * 100:.3f}% of the {error["limit"] * 1e3:.0f} mHa cap
    and <code>E_gap</code> clears its strict floor by
    {e_gap_slack:.2f}&nbsp;×&nbsp;10⁻⁶ Ha — {len(pinned)} of the
    {len(margins)} constraints spent to within {(1 - PINNED_FRACTION) * 100:.1f}% of their
    limit — and the widest ANGLE circuit is {int(angle["value"])} of {int(angle["limit"])}
    qubits, {angle["fraction_of_limit_used"] * 100:.1f}% of its cap. The
    {self.plateau_hi - self.median:.1f}-wide band is what that leaves.</p>
    <figure>
      {headroom_svg(margins)}
      <figcaption>Each hard constraint on the <code>dfthcsos</code> task, as a fraction of its
      limit, measured by our own verifier over the submitted files and recorded in
      <code>{PARTS_PATH.name}</code>. Amber is pinned, meaning within
      {(1 - PINNED_FRACTION) * 100:.1f}% of the limit: {" and ".join(pinned)}. The
      ANGLE width has {int(angle["limit"]) - int(angle["value"])} qubits of room. The one-norm
      is plotted against the floor that binds it, {min_lcu_lambda(FEMOCO_C0):.2f}, and not
      against <code>MIN_LAMBDA</code>: that rule floors Λ rather than λ_eff — this
      submission&rsquo;s λ_eff of {plateau["lambda_eff"]:.3f} is below
      {lcu["limit"]:g} and still passes — but at {lcu["limit"]:g} it never binds, and plotting
      it there drew Λ with 84% of room it does not have.</figcaption>
    </figure>
    <p>The move this part made is now legible. Its shape is <strong>R={r}, B={b}, C={c}</strong>
    — its author&rsquo;s report, and the shape of the arrays in the <code>submission.npz</code>
    we submitted — against the baseline&rsquo;s <strong>R={br}, B={bb}, C={bc}</strong>, which
    the verifier declares itself. So the plateau <strong>cut the rank from {br} to {r} and
    widened the basis from {bb} to {b}</strong>: {ref_vectors:,}&nbsp;→&nbsp;{vectors:,} LCU
    vectors for the ANGLE circuits to carry (the verifier builds <code>N + R·B</code> of them),
    Toffolis {reference.toffolis:,.0f}&nbsp;→&nbsp;{plateau["T_total"]:,.0f}, width
    {int(reference.detail["num_qubits"])}&nbsp;→&nbsp;{plateau["Q_total"]}. What it spent is
    ensemble capacity, R·C, down {ref_capacity / capacity:.1f}× from {ref_capacity:,} to
    {capacity:,} — and the chart above is where that shows: the accuracy cap the reference
    clears at {ref_error * 100:.0f}% is the one this part is sitting on.
    <strong>Low rank is what the plateau found, and the accuracy cap is what stops it going
    lower.</strong> Searching for a cheaper shape than R·B&nbsp;=&nbsp;{r * b:,} was to be the
    next experiment. It is closed: λ_eff cannot fall below
    {lambda_eff_floor(FEMOCO_C0):.3f} for <em>any</em> factorization, against this
    submission&rsquo;s {plateau["lambda_eff"]:.3f} — so the whole search is worth
    {(1 - lambda_eff_floor(FEMOCO_C0) / plateau["lambda_eff"]) * 100:.2f}% of one of three cost
    factors, about +0.24 of the +1.82 the blend needs. Why, and what is left, is below.</p>
"""

    def what_the_page_omits(self) -> str:
        invsqrt, dqi = self.reports["invsqrt"], self.reports["dqi"]
        rotation_eps = float(invsqrt.detail["rotation_eps"])
        return f"""
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
          <td class="num">{int(invsqrt.detail["rotation_count"])} rotations at
          ε={rotation_eps:.2e}<br>{rotation_t_count(rotation_eps):,.2f} T gates each</td></tr>
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
"""

    def repetition_floor(self) -> str:
        dqi = self.reports["dqi"]
        overhead = (dqi.repetition_factor - 1) * 100
        return f"""
    <h2>The {overhead:.0f}% nobody can optimise away</h2>
    <p class="finding">A decoder that fails zero of {dqi.shots:,} shots is still bounded at
    ε₊ = {dqi.eps_plus:.2e}, so <code>dqi</code> pays R = {dqi.repetition_factor:.4f} — a
    {overhead:.0f}% repetition overhead — before any circuit work. Only more shots lower that floor, and the shot
    count is the verifier&rsquo;s to choose.</p>
    <figure>
      {repetition_svg(dqi)}
      <figcaption>R(ε) = ln2 / −ln(½ + 2√ε), drawn from
      <code>cost.dqi_repetition_factor</code>. The cost diverges at ε = 1/16, where the verifier
      rejects outright. Our point is the Clopper–Pearson bound for 0 failures in
      {dqi.shots:,} shots.</figcaption>
    </figure>
"""

    def what_chemistry_controls(self) -> str:
        dfthcsos = self.reports["dfthcsos"]
        resources = dfthcsos.detail["resources"]
        # The reference's own headroom, off the same constraint table as the plateau part's.
        reference = dfthcsos_constraint_margins(dfthcsos.detail)
        angle, error, e_gap = (
            reference["angle_qubits"],
            reference["hamiltonian_error_ha"],
            reference["e_gap_ha"],
        )
        reference_res = ShapeResources.from_report(dfthcsos.detail)
        plateau = self.parts["parts"]["dfthcsos"]
        plateau_free = float(plateau["T_angles"])
        plateau_bought = float(plateau["T_total"]) - plateau_free
        blend = self.our_score
        cut = cost_cut_for_blend(blend, ESCAPE_SCORE)
        return f"""
    <h2>What a chemistry submission actually controls</h2>
    <p class="finding">The ANGLE circuits a submission ships are
    {resources["angles_toffoli"] / resources["toffoli"] * 100:.0f}% of the Toffoli count but only
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
    <h2>What one Toffoli is worth, and which ones are for sale</h2>
    <p class="finding">Driving the verifier&rsquo;s own cost model at the plateau&rsquo;s shape
    reproduces its recorded T and Q exactly, and splits them: <strong>{plateau_free:,.0f} of the
    {plateau["T_total"]:,.0f} Toffolis ({plateau_free / plateau["T_total"]:.1%}) are ANGLE
    circuits</strong>, which a submission writes and so can make cheaper at <em>identical</em>
    accuracy. The rest is closed-form in <code>(R, B, C)</code> and <code>brot</code>, payable
    only in accuracy &mdash; and the accuracy is spent, the binding residual at
    {plateau["constraint_margins"]["hamiltonian_error_ha"]["fraction_of_limit_used"] * 100:.3f}%
    of its cap.</p>
    <figure>
      {toffoli_split_svg([
          ("reference", reference_res.angle_toffoli, float(reference_res.heuristic_toffoli)),
          ("plateau (submitted)", plateau_free, plateau_bought),
      ])}
      <figcaption>The plateau cut the total 5.7&times; and, in doing so, cut the free share from
      {reference_res.angle_share:.0%} to {plateau_free / plateau["T_total"]:.0%} &mdash; the
      heuristic half shrank faster than the loaders did. Reference split read from its own
      verifier report; the plateau&rsquo;s from
      <code>T_angles</code>&nbsp;/&nbsp;<code>T_total</code> in its parts record.</figcaption>
    </figure>
    <p>At this cost, one Toffoli is worth
    <strong>{marginal_blend_per_unit(blend, plateau["T_total"]):+.6f}</strong> blend and
    one qubit {marginal_blend_per_unit(blend, plateau["Q_total"]):+.6f} &mdash; a qubit
    is worth {plateau["T_total"] / plateau["Q_total"]:.0f} Toffolis purely because Q is that much
    smaller. Reaching {ESCAPE_SCORE:g} needs
    <strong>{cut:.2%}</strong> off this cost:
    {plateau["T_total"] * cut:,.0f} Toffolis, or
    {plateau["Q_total"] * cut:,.0f} qubits, or any
    mix. <strong>Out of the free share alone that is loaders
    {plateau["T_total"] * cut / plateau_free:.0%}
    cheaper</strong> &mdash; large, and the only direction left that costs no accuracy. What
    those circuits turn out to be is what changes that number&rsquo;s meaning.</p>
    <p class="note">Every way of buying the bought half is closed, measured rather than argued.
    <code>brot</code> cannot go down: re-rounding the submitted angles onto the
    <code>brot&nbsp;=&nbsp;7</code> grid breaches the 20&nbsp;mHa cap by <strong>118&times;</strong>
    on <code>error_1</code>, and they already sit on that cap at 8, so a fresh fit at half the
    grid spacing would have to beat re-rounding by two orders of magnitude.
    <code>brot</code> cannot go up: at 9 the inner preparation&rsquo;s QROAM jumps to a branch
    taking the shared workspace from 304 to 607 qubits, costing <strong>&minus;3.27 blend</strong>
    despite 1,316 fewer Toffolis. And the shape cannot be re-split: at a constant
    R&middot;B&nbsp;=&nbsp;1,144 the submitted (8,&nbsp;143) beats all 16 integer factor pairs.</p>
    <p class="note">The reference sits close to three cliffs, all of which reject outright rather
    than gracefully: the widest ANGLE circuit at {int(angle["value"])} qubits against the
    {int(angle["limit"])} cap ({angle["fraction_of_limit_used"] * 100:.0f}% of it), the binding
    accuracy term at {error["value"]:.5f} Ha against a {error["limit"]:.3f} Ha cap
    ({error["fraction_of_limit_used"] * 100:.0f}%), and <code>E_gap</code> at
    {e_gap["value"]:.4f} Ha against a strict {e_gap["limit"]} Ha floor
    (+{(e_gap["value"] / e_gap["limit"] - 1) * 100:.1f}%). Those are the
    <em>reference&rsquo;s</em> margins, read off the same constraint table; the chart above
    measures the same three for the plateau part, which has spent nearly all of them.</p>
"""

    def loader_floor(self) -> str:
        plateau = self.parts["parts"]["dfthcsos"]
        addresses = lcu_vectors(int(plateau["R"]), int(plateau["B"]))
        resources = ShapeResources.from_parts_record(plateau)
        sweep = AngleFloorSweep.over(
            resources, int(plateau["brot"]), addresses, self.our_score
        )
        submitted, per_pass = sweep.submitted, resources.angle_toffoli / 2
        best_floor, best_rate = sweep.best_floor, sweep.best_rate
        target = sweep.at(best_rate.batch_size)
        paired = best_floor.at_floor * blend_factor(
            lambda_eff_floor(FEMOCO_C0), float(plateau["lambda_eff"])
        )
        return f"""
    <h2>The loaders are Boolean circuits, and their floor still clears the band</h2>
    <p class="finding">The free half is not a lookup table. Each <code>ANGLE_k</code> is a GF(2)
    circuit: a ladder of AND-ancillas, {LOADER_REBUILDS["circuit_cx"]:,} CNOTs fanning their
    products into the word bits &mdash;
    <strong>{LOADER_REBUILDS["selection_cx"]:,} of them targeting the selection register
    itself</strong>, used in place as a linear scratchpad and restored by the end &mdash; then
    measurement and a conditioned <code>CZ</code> clearing every ancilla. The verifier charges
    Toffoli-like gates only, and a conditioned gate only on the shots whose bit is set, so
    <strong>XOR costs nothing and the AND gates are the whole bill</strong>. A loader&rsquo;s
    price is the multiplicative complexity of its word table, not its address count &mdash;
    which is why {addresses:,} addresses cost
    {LOADER_REBUILDS["cheapest_circuit_toffoli"]} Toffolis and not {addresses - 1:,}.</p>
    <p>Rebuilding the same <code>angle_words</code> three ways sizes the construction:
    kickmix&rsquo;s own <code>init_lookup</code> costs
    <strong>{LOADER_REBUILDS["init_lookup_toffoli"]:,}</strong> Toffolis per pass over the
    family, the starter kit&rsquo;s segment-tree QROM
    <strong>{LOADER_REBUILDS["qrom_toffoli"]:,}</strong>, and the submitted family
    <strong>{per_pass:,.0f}</strong> &mdash;
    <strong>{LOADER_REBUILDS["init_lookup_toffoli"] / per_pass:.1f}&times; cheaper than the
    library primitive</strong> and {LOADER_REBUILDS["qrom_toffoli"] / per_pass:.1f}&times;
    cheaper than the starter kit, bought with
    {LOADER_REBUILDS["extra_workspace_qubits"]} extra workspace qubits. It is nobody&rsquo;s
    default.</p>
    <figure>
      {angle_floor_svg(sweep)}
      <figcaption>A straight-line GF(2) circuit with <em>n</em> inputs, <em>A</em> ANDs and
      <em>b</em> outputs realises at most 2<sup>2An+A&sup2;+A+b(n+A+1)</sup> maps, and must hit a
      given <em>b</em>-bit word on each of {addresses:,} addresses, so almost every table needs
      A&nbsp;&ge;&nbsp;{submitted.and_floor} here. <em>n</em> is {submitted.public_qubits}, not
      {submitted.address_bits}: each circuit starts with its word registers holding the previous
      batch, itself a function of the address, and the submitted circuits use them. Withholding
      them would give {submitted.address_only_floor} and a floor that never clears the band
      &mdash; a bound must allow every resource a construction may legitimately use.</figcaption>
    </figure>
    <p>The submitted mean of {sweep.measured_ands:,.1f} ANDs per circuit is
    <strong>{sweep.rate:.3f}&times;</strong> that floor, and <code>J</code> &mdash; the batch
    size &mdash; is a free parameter: the contract fixes only that a family covers all
    {NUM_ANGLE_WORDS} angles and that no circuit exceeds {MAX_ANGLE_QUBITS} qubits, so
    <em>K</em> falls as 1/<em>J</em> while a circuit&rsquo;s floor grows only as &radic;<em>J</em>.
    Two numbers follow. <strong>Re-batching alone is
    {best_rate.at_rate - sweep.blend:+.3f} blend</strong>, to {best_rate.at_rate:.5f} at
    J={best_rate.batch_size}, with no new circuit cleverness and no accuracy change. And
    reaching {sweep.target:g} needs <strong>{target.ands_for_target:,.0f} ANDs per circuit at
    J={target.batch_size}, {target.target_over_floor:.2f}&times; the floor</strong> &mdash;
    where this family already builds {LOADER_REBUILDS["cheapest_circuit_toffoli"]}&ndash;{
        LOADER_REBUILDS["dearest_circuit_toffoli"]} per circuit at J={submitted.batch_size}.
    The target is a J={target.batch_size} construction no dearer per circuit than the ones
    already written, over
    {target.family.output_bits / submitted.output_bits:.1f}&times; the output bits.</p>
    <p class="note">So the loaders are the one direction not closed: the floor reaches
    {best_floor.at_floor:.3f} at J={best_floor.batch_size}, above the band, and with
    <code>lambda_eff</code> at its own floor too the pair reach {paired:.2f}. What is not
    established is that any such construction exists. Going from {sweep.rate:.2f}&times; to
    {target.target_over_floor:.2f}&times; is Boolean-minimisation work this project has not
    done, the {target.family.free_workspace} workspace qubits a J={target.batch_size} batching
    leaves are untested, and a counting argument bounds almost every table rather than proving
    anything about this one &mdash; it says what is possible, not what is reachable.</p>
"""

    def cost_anatomy(self) -> str:
        invsqrt, dqi = self.reports["invsqrt"], self.reports["dqi"]
        toffoli_span = dqi.toffolis / invsqrt.toffolis
        return f"""
    <h2>Cost anatomy</h2>
    <p class="finding">A qubit saved is worth twice as much on <code>dqi</code>, whose cost carries
    Q²; the other two are linear in Q. That exponent, not the raw counts, decides which lever
    pays — and the three tasks span {toffoli_span:,.0f}× in Toffoli count, so they are not
    interchangeable targets.</p>
    <div class="tw"><table>
      <thead><tr><th>Task</th><th>Reference cost</th><th>Scored as</th><th>Q</th>
      <th>Toffolis/shot</th><th>Leading factor</th><th>Shots</th></tr></thead>
      <tbody>{self.reference_rows()}</tbody>
    </table></div>
    <p class="note">Each row&rsquo;s cost is reconstructable from the columns beside it. For
    <code>invsqrt</code> the leading factor
    C = {float(invsqrt.detail["toffoli_equivalent"]):,.2f}
    already contains the Toffolis, so its cost is Q × C and the Toffoli column is shown for
    comparison only.</p>
"""

    def local_loop_cost(self) -> str:
        return f"""
    <h2>What the local loop costs</h2>
    <p class="finding">A full three-task verification takes about
    {sum(HAND_RECORDED["minutes"].values()):.0f} minutes on 10 CPU cores and is free and
    unlimited — against {HAND_RECORDED["submissions_per_day"]} rationed submissions a day. That
    ratio is the whole reason the dev set is worth building.</p>
    <div class="funnel">{self.run_log_rows()}</div>
    <p class="note">Not a funnel yet: nothing has been eliminated, because no candidates have been
    generated. The per-stage denominators this section will carry only exist once a generator is
    running.</p>
"""

    def seed_noise(self) -> str:
        spread, seeds = self.seed_spread, self.seed_runs["seeds"]
        shots = self.seed_runs["shots_per_run"]
        clean = next(s for s, f in zip(seeds, self.failed_shots) if not f)
        failed = next(s for s, f in zip(seeds, self.failed_shots) if f)
        hits = sum(1 for f in self.failed_shots if f)
        lo, hi = self.failure_interval
        rate = sum(self.failed_shots) / (len(seeds) * shots)
        return f"""
    <h2>How much of a difference is the seed?</h2>
    <p class="finding">Re-verifying one unchanged submission under {len(seeds)} verifier seeds
    splits the three tasks into two regimes. <strong>Chemistry is deterministic and
    <code>dqi</code> is not</strong>: <code>dfthcsos</code> returned the identical cost every
    seed and <code>invsqrt</code> moved by {spread["invsqrt"]["ppm"]:.2f} ppm, while
    <code>dqi</code> moved by {spread["dqi"]["ppm"]:,.0f} ppm and carries the blend&rsquo;s
    whole {spread["blend"]["range"]:.4f}-point spread.</p>
    <div class="tw"><table>
      <thead><tr><th></th><th class="num">min</th><th class="num">max</th>
      <th class="num">spread</th><th>reading</th></tr></thead>
      <tbody>{self.seed_rows()}</tbody>
    </table></div>
    <p>The mechanism is one discrete step, not a distribution.
    {len(seeds) - hits} of {len(seeds)} seeds recorded
    <code>empirical_error_rate = 0</code> and returned identical numbers to every digit &mdash;
    a one-sided bound of {clean["dqi_detail"]["effective_error_rate"]:.4e} and a repetition
    factor of {clean["dqi_detail"]["repetition_factor"]:.5f}. Seed {failed["seed"]} recorded
    <code>{failed["dqi_detail"]["empirical_error_rate"]:.6e}</code>, exactly
    1&nbsp;/&nbsp;{shots:,} &mdash; <strong>one failed shot</strong> &mdash; lifting the bound
    to {failed["dqi_detail"]["effective_error_rate"]:.4e} and the factor to
    {failed["dqi_detail"]["repetition_factor"]:.5f}, a
    {(failed["dqi_detail"]["repetition_factor"] / clean["dqi_detail"]["repetition_factor"] - 1) * 100:.3f}%
    cost penalty. Nothing else moved with it: the simulated per-shot Toffoli mean varied by
    {spread["toffolis"]["ppm"]:.2f} ppm and <code>num_qubits</code> held at
    {clean["dqi_detail"]["num_qubits"]:,}.</p>
    <p class="note">So {len(seeds)} seeds size that step precisely and its frequency badly.
    {sum(self.failed_shots)} failure in {len(seeds) * shots:,} shots puts the per-shot rate near
    {rate:.2e}, which would fail a {shots:,}-shot run about
    {1 - math.exp(-rate * shots):.0%} of the time &mdash; but the 95% interval on that
    frequency, from {hits} of {len(seeds)} runs, is [{lo:.1%}, {hi:.1%}]. The leaderboard draws
    its own seed, so our {self.our_score:.5f} was a zero-failure draw and the identical file
    could legitimately return {spread["blend"]["min"]:.4f}.</p>
"""

    def blend_check(self) -> str:
        public = score(PUBLIC_COSTS, DEFAULT_CALIBRATION)
        delta_ppm = abs(public - PUBLIC_BLEND) / PUBLIC_BLEND * 1e6
        floor_ppm = blend_precision_ppm(PUBLIC_COSTS, PUBLIC_COST_HALF_ULP)
        return f"""
    <h2>Checking the blend arithmetic</h2>
    <p class="finding">Feeding the public frontier&rsquo;s self-reported per-task costs through our
    scoring path reproduces the blend they reported, {PUBLIC_BLEND}, to
    {delta_ppm:.1f} parts per million. That checks the ratio-and-geometric-mean
    arithmetic against a third party&rsquo;s independent calculation — it does <em>not</em> check
    the three per-task cost formulas, since both sides use the same self-reported costs.</p>
    <div class="tw"><table>
      <thead><tr><th>Task</th><th>Public best cost (third-party)</th>
      <th>Ratio vs our reference</th></tr></thead>
      <tbody>{self.ratio_rows()}</tbody>
    </table></div>
    <p class="note">Computed with <code>cost.task_ratio</code> and <code>cost.score</code> against
    the measured calibration, giving {public:.5f}. The agreement cannot be read as
    finer than the inputs: the rounding in those costs as quoted moves the blend by about
    ±{floor_ppm:.1f} ppm on its own, so {delta_ppm:.1f} ppm means &ldquo;agrees to the
    precision available&rdquo;, not a measured resolution.</p>
"""

    def verdict_ledger(self) -> str:
        return f"""
    <h2>Verdict ledger</h2>
    {self.verdict_blocks()}
    <h3>Promotion gates</h3>
    <div class="tw"><table><tbody>{self.gate_rows()}</tbody></table></div>
"""

    def not_established(self) -> str:
        plateau = self.parts["parts"]["dfthcsos"]
        return f"""
    <h2>Not established</h2>
    <ul>
      <li><strong>We have not compressed anything ourselves.</strong> The
      {self.our_score:.5f} is an assembly of other people&rsquo;s published circuits; everything we
      built is measurement apparatus. No claim here is a claim about a circuit we designed.</li>
      <li><strong>How often <code>dqi</code> drops a shot is unmeasured.</strong> Its
      size is now measured — one failed shot, worth
      {self.seed_spread["blend"]["range"]:.4f} blend points — but
      {sum(self.failed_shots)} failure across {len(self.failed_shots)} seeds leaves the
      <em>frequency</em> anywhere in
      [{self.failure_interval[0]:.1%}, {self.failure_interval[1]:.1%}], so the expected penalty
      on any one submission is not known. Pinning it needs the &ge;&nbsp;48 draws gate G2
      asks for.</li>
      <li><strong>Every third-party number here is unverified by us.</strong> The per-task
      frontier costs and the {PUBLIC_BLEND} blend are self-reported in Kaggle topic 754174, and
      the blend check above rests entirely on them — if they are wrong, that check says nothing.
      The leaderboard scores and the {self.standing.teams:,} team count are the competition&rsquo;s
      own. So are the <code>R</code>/<code>B</code>/<code>C</code> and λ_eff figures in
      <code>{PARTS_PATH.name}</code> — those are each part author&rsquo;s own report, copied from
      the dataset manifest; only <code>our_measurement</code>, the constraint margins and
      <code>J</code> are ours.</li>
      <li><strong>Whether <code>T</code> and <code>Q</code> have any headroom is
      unmeasured.</strong> The λ_eff floor says where the chemistry score is <em>not</em> — no
      factorization, of any shape, moves that factor by more than
      {(1 - lambda_eff_floor(FEMOCO_C0) / plateau["lambda_eff"]) * 100:.2f}% — but it says
      nothing about the loaders, whose own floor is a counting bound on almost every table
      rather than a statement about this one. Whether a construction reaching it exists is the
      open measurement, and it may come back out of reach.</li>
      <li><strong>Two leaderboard readings, neither on a circuit we designed.</strong> The loop
      agreed with the leaderboard on the reference and on this assembly of published circuits.
      Agreement there does not prove agreement on a candidate our own search produces.</li>
      <li><strong>Which task to commit to is open.</strong> The structural facts above bear on it;
      the decision does not follow from them.</li>
      <li><strong>Why {self.escaped} teams cleared {ESCAPE_SCORE:g} while half the field sits near
      {self.median:.1f}</strong> is unexplained, and which task holds that margin is the most
      valuable unknown in the competition.</li>
    </ul>
"""

    def footer(self) -> str:
        return f"""
  <footer>
    Built by <code>scripts/build_artifact.py</code> from the three official verifier outputs
    committed under <code>tests/data/</code>, the verdict ledger and leaderboard snapshot in
    <code>docs/experiments/</code>, and the formulas and constants imported from
    <code>quantum_compilation.cost</code>. The three wall-clock run times and the
    {HAND_RECORDED["submissions_per_day"]}-a-day submission budget are recorded by hand and live
    in <code>HAND_RECORDED</code>, the ANGLE circuit readings and loader rebuilds in
    <code>LOADER_REBUILDS</code>; every ratio between them is computed, and third-party figures
    are labelled where they appear. Status
    {STATUS_DATE}.
  </footer>
"""

    def html(self) -> str:
        # Everything between the header and the footer is a <section>, so the tag is written
        # once here rather than at both ends of all thirteen.
        body = "".join(
            f"\n  <section>{section()}  </section>\n"
            for section in (
                self.strip, self.where_we_stand, self.why_the_plateau_is_a_wall,
                self.what_the_page_omits, self.repetition_floor, self.what_chemistry_controls,
                self.loader_floor, self.cost_anatomy, self.local_loop_cost,
                self.seed_noise, self.blend_check, self.verdict_ledger,
                self.not_established,
            )
        )
        sections = f"{self.header()}{body}{self.footer()}"
        return f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>VQCC Compression Log</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>
{STYLE}</style>

<main>{sections}</main>
</html>
"""


def main() -> None:
    artifact = Artifact.load()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(artifact.html())
    print(f"wrote {OUTPUT_PATH} ({OUTPUT_PATH.stat().st_size:,} bytes)")
    print(f"wrote {write_progress_png(artifact.verdicts)}")


if __name__ == "__main__":
    main()
