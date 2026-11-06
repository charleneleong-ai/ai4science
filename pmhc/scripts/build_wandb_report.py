"""Build and refresh this project's one W&B report.

The narrative sections — Highlights, Task, Approach, What we learned, Next — are written here
once and survive every refresh. Results and the full run log are regenerated in place from
`verdicts.yaml` and from the runs `log_wandb.py` wrote, so a refresh never overwrites the prose
and the prose never goes stale against the ledger. Every number in both halves is read from
`campaign.yaml`, so nothing on the report is typed in twice.

    uv run --extra wandb python scripts/build_wandb_report.py create --entity chaleong
    uv run --extra wandb python scripts/build_wandb_report.py refresh --url <report url>
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import typer
from build_artifact import grid_collapse, read_yaml

app = typer.Typer(add_completion=False)

TITLE = "Reading a neoepitope — KRAS G12D on HLA-A*11:01"
BRANCH = "feat/pmhc-binder"
RESULTS = "Results"
RUN_LOG = "Full run log"
REPO = "https://github.com/charleneleong-ai/ai4science"
PR = f"{REPO}/pull"
BLOB = f"{REPO}/blob/{BRANCH}/pmhc"


def api_overrides() -> dict[str, str]:
    """Never assume wandb.ai: the deployment may be self-hosted."""
    host = os.environ.get("WANDB_BASE_URL")
    return {"base_url": host} if host else {}


def novelty_removed(data: dict[str, Any]) -> int:
    """What the two novelty gates took out, before anything was scored for binding."""
    return sum(s["removed"] for s in data["funnel"] if s["stage"].endswith("novelty"))


def highlights(data: dict[str, Any]) -> str:
    head = data["headline"]
    cone = data["escape_cone"]
    var = data["variance"]["designs"]
    return f"""## Highlights

**{data["updated"]} — {head["submitted"]} designs selected for the 20 slots; not yet uploaded, and
no external score.** Adaptyv scores after the {data["deadline_aoe"]} AoE close, and the portal
takes entries by hand, so every number here is this project's own measurement against
its own calibrated bar, and the bar is a real molecule rather than a threshold.

| what | value |
|---|---|
| best design's peptide-dependence drop | `{head["best_design"]}` {float(head["best_drop"]):+.3f} over runs {head["best_per_run"]} |
| the bar: the 7OW6 receptor's own drop | {float(head["bar"]):+.3f} over {data["bar"]["invocations"]} byte-identical invocations |
| designs at or above the bar | {head["above_bar"]} above, {head["tied_bar"]} tied, {head["below_bar"]} below, of {head["screened"]} screened |
| designs clearing zero at all | {head["clear_zero"]} of {head["screened"]} |
| **ranking tier 1, on the control** | **{float(head["tier1_control_bias"]):+.4f} ± {float(head["tier1_control_bias_sem"]):.4f}** over {head["tier1_control_invocations"]} invocations — the wrong sign |

**The levers that moved the result, and the one that did not.**

1. **The mutation is chemically unreachable, so the strategy is indirect.** The G12D carboxylate
   oxygen has a {cone["by_reach_angstrom"][8]}% escape cone at 8 Å reach
   ({cone["free_rays_at_8A"]} of {cone["rays"]:,} rays), collapsing
   {grid_collapse(data)}-fold between 6 and 7 Å. Designs target the conformational bulge at peptide positions 5–8 instead.
2. **The denominator was the lever, not the generator.** Two runs per design cannot estimate that
   design's own between-run SD, so verdicts divide by the cohort-pooled one: between-run
   {var["between_sd"]} against within-run {var["within_sem"]}, a run effect of
   {var["run_effect"]}. On that denominator {head["above_bar"]} designs clear the bar and
   {head["tied_bar"]} tie it — a count a two-sample error bar was never able to support.
3. **Novelty had to gate before score, not after.** ipTM rewards handing the target its own chains
   back, so a score-first ranking would have promoted plagiarised designs. The two novelty gates
   removed {novelty_removed(data)} of {data["funnel"][0]["surviving"]} before any binding
   score was read.
4. **Ranking tier 1 never worked, and is reported rather than gated.** Measured on the receptor
   that is specific by construction, the predictor prefers the peptide it does *not* recognise.
   Selectivity in this submission is sampled, never engineered, and the methodology says so.

**Progression — the funnel, with every denominator.**

| stage | surviving | removed | why |
|---|---|---|---|
""" + "\n".join(
        f"| {s['stage']} | {s['surviving']} | {s['removed']} | {s['reason']} |"
        for s in data["funnel"]
    )


def task(data: dict[str, Any]) -> str:
    ind = data["indirect_readout"]
    return f"""## Task

Design a binder for the KRAS G12D neoepitope `VVVGADGVGK` presented on HLA-A\\*11:01 that
prefers it over the wild-type `VVVGAGGVGK` — one aspartate against one glycine, at peptide
position 6.

**The metric.** A design's *peptide-dependence drop* is its complex score minus the same
complex scored with the peptide removed, so a binder that reads the peptide drops and one that
binds the HLA groove does not. The bar is the 7OW6 receptor's own drop, {float(data["headline"]["bar"]):+.3f}.

**Why selectivity cannot be gated here.** The aspartate points into the groove, so the two
targets differ at the binder's surface only by a backbone shift of
{ind["backbone_shift_angstrom"][6]}/{ind["backbone_shift_angstrom"][7]}/{ind["backbone_shift_angstrom"][8]} Å
at positions 6/7/8, against a {ind["superposition_residual"]} Å superposition residual — at or
below the predictor's own coordinate error.
"""


def approach() -> str:
    return f"""## Approach

The pipeline as it stands, in the order a design meets it:

1. **Generate** against the intact 7OW6 complex with the TCR deleted, 64 designs for each of four
   hotspot variants (`bulge_core`, `bulge_wide`, `mut_proximal`, `hub`), all aimed at peptide
   positions 5–8.
2. **Sequence novelty** by local alignment, not global: the global metric scored a verbatim
   40-residue slice of the HLA heavy chain at 0.000 and would have passed it as novel.
3. **Fold novelty** by FoldSeek TM-score against the target's own chains and the PDB.
4. **Engagement**: the peptide's share of the design's interface must reach the receptor's own
   0.200, and the nominated hotspot must actually be touched.
5. **Score** the peptide-dependence drop, replicated across independent Boltz invocations, and
   rank against the cohort-pooled between-run SD.
6. **Report** the wild-type difference for every finalist without gating on it.

Entry points: [`scripts/cofold_followup.py`]({BLOB}/scripts/cofold_followup.py)
for the statistics, [`scripts/submission.py`]({BLOB}/scripts/submission.py)
for the slot assignment, [`scripts/build_artifact.py`]({BLOB}/scripts/build_artifact.py)
for the page. The bar's calibration is in the ledger below as `bar-reproducibility`.
"""


def learned(data: dict[str, Any]) -> str:
    bar = data["bar"]
    return f"""## What we learned

**Worked, with the measured gain.**

| lever | measured effect |
|---|---|
| local-alignment novelty metric | caught the {data["funnel"][1]["removed"]} verbatim β2-microglobulin copies the global metric scored 0.000 |
| novelty before score | {novelty_removed(data)} of {data["funnel"][0]["surviving"]} removed before any binding score was read |
| cohort-pooled denominator | a two-sample SD carries ~70% relative error, so every verdict divides by the pooled between-run SD instead |
| replication of the bar | {bar["invocations"]} invocations put it at {float(bar["pooled"]):+.3f} ± {bar["pooled_spread"]}, against {float(bar["cohort_run1"]):+.3f} from the single run that was later excluded on structural grounds |
| pairing both peptides inside one invocation | the run offset cancels exactly; taken from different invocations it does not, and is not in the error bar either |

**Closed — do not retry.**

- *Gating on the G12D atom itself.* No approach vector exists; the escape cone is
  {data["escape_cone"]["by_reach_angstrom"][8]}% at 8 Å.
- *Ranking tier 1 as a selectivity gate.* Inverted on the control at
  {float(data["headline"]["tier1_control_bias"]):+.4f}; it measures the predictor, not the design.
- *Per-design between-run SD from two runs.* ~70% relative error; it is not an error bar.
- *ipTM as a first-pass ranking.* It prefers the target's own chains to any de novo design.
"""


def nxt(data: dict[str, Any]) -> str:
    return f"""## Next

- **Pending:** the upload itself. The 20-slot submission is assembled and validated, but the
  Proteinbase portal takes entries by hand and has no API, so it is not in until it is uploaded.
- **Queued, and not met:** the second-predictor pass. The standing discipline on this project is
  that every load-bearing claim gets an independent method, and the Chai-1 container was never
  located, so every number here is Boltz 2.2.1 scoring structures Boltz 2.2.1 predicted.
  Replication bounds sampling and run-to-run error; it cannot detect an error shared by every
  run, and the inverted tier-1 result is exactly that kind of error.
- **Deadline:** {data["deadline_aoe"]} AoE, {data["deadline_sydney"]} Sydney.
"""


def entries(ledger: dict[str, Any]) -> list[dict[str, Any]]:
    return [{**e, "label": e["id"], "result": str(e["result"]).strip()} for e in ledger["experiments"]]


def results_markdown(data: dict[str, Any], ledger: dict[str, Any]) -> str:
    from autoresearch.wandb_report import verdict_table  # only in the `wandb` extra

    head = data["headline"]
    table = verdict_table(entries(ledger), pr_url=PR)
    return f"""## {RESULTS}

**External scores: none yet.** Adaptyv scores after the {data["deadline_aoe"]} AoE close. Until
then the only measurements are local, against the receptor's own bar of {float(head["bar"]):+.3f}.

**Verdict ledger** — one row per experiment, its gate, and what it settled. Gates are defined
once at the top of
[`docs/experiments/pmhc/verdicts.yaml`]({BLOB}/docs/experiments/pmhc/verdicts.yaml).

{table}
The panels below are the committed measurement tables, carried by the run each belongs to.
"""


def run_log_markdown(ledger: dict[str, Any]) -> str:
    rows = "\n".join(
        f"| `{e['id']}` | {e['verdict']} | {str(e['what']).strip()} | {str(e.get('consequence', '')).strip()} |"
        for e in ledger["experiments"]
    )
    not_established = "\n".join(f"- {str(item).strip()}" for item in ledger["not_established"])
    return (
        f"## {RUN_LOG}\n\n| experiment | verdict | what it tested | consequence |\n|---|---|---|---|\n"
        f"{rows}\n\n**What is not established**\n\n{not_established}\n"
    )


def table_panels(entity: str, project: str) -> list[Any]:
    """One panel per `wandb.Table` actually logged, read off the runs rather than hardcoded, so a
    key renamed at log time shows up as a missing panel instead of a silently empty one."""
    import wandb
    import wandb_workspaces.reports.v2 as wr

    api = wandb.Api(overrides=api_overrides())
    # lazy=False fetches summaryMetrics inline: the default issues the lightweight fragment,
    # which omits them, so every `run.summary` below would be its own round trip.
    runs = api.runs(f"{entity}/{project}", filters={"tags": "pmhc"}, lazy=False)
    keys = sorted(
        {key for run in runs for key in run.summary.keys() if key.startswith("table/")}
    )
    return [wr.WeavePanelSummaryTable(table_name=key) for key in keys]


@app.command()
def create(
    entity: str = typer.Option(..., help="W&B entity"),
    project: str = typer.Option("pmhc-g12d-binder"),
    campaign: Path = typer.Option(Path("docs/experiments/pmhc/campaign.yaml")),
    verdicts: Path = typer.Option(Path("docs/experiments/pmhc/verdicts.yaml")),
    artifact_url: str = typer.Option("", help="the published Claude Artifact, linked at the top"),
) -> None:
    """Write the report once, narrative and all. Refreshes go through `refresh`.

    Panels are bound to whatever tables the runs already carry, so this is run after
    `log_wandb.py`; run before it, the report lands without them and a `refresh` adds them.
    """
    import wandb_workspaces.reports.v2 as wr

    data = read_yaml(campaign)
    ledger = read_yaml(verdicts)
    link = f" · [the illustrated write-up]({artifact_url})" if artifact_url else ""
    blocks: list[Any] = [
        wr.TableOfContents(),
        wr.CalloutBlock(
            text=f"{data['task']}. Code and full commit history in [PR #103]({PR}/103){link}."
        ),
    ]
    blocks += [
        wr.MarkdownBlock(text=text)
        for text in (
            highlights(data),
            task(data),
            approach(),
            results_markdown(data, ledger),
        )
    ]
    panels = table_panels(entity, project)
    if panels:
        blocks.append(
            wr.PanelGrid(runsets=[wr.Runset(entity=entity, project=project)], panels=panels)
        )
    blocks.append(
        wr.H2(text=RUN_LOG, collapsed_blocks=[wr.MarkdownBlock(text=run_log_markdown(ledger))])
    )
    blocks += [wr.MarkdownBlock(text=text) for text in (learned(data), nxt(data))]

    report = wr.Report(
        project=project, entity=entity, title=TITLE, description=data["task"], blocks=blocks
    )
    report.save()
    typer.echo(report.url)


@app.command()
def refresh(
    url: str = typer.Option(..., help="the report's URL"),
    entity: str = typer.Option(..., help="W&B entity"),
    project: str = typer.Option("pmhc-g12d-binder"),
    campaign: Path = typer.Option(Path("docs/experiments/pmhc/campaign.yaml")),
    verdicts: Path = typer.Option(Path("docs/experiments/pmhc/verdicts.yaml")),
) -> None:
    """Rewrite Results and the collapsed run log in place; the narrative sections are untouched."""
    from autoresearch.wandb_report import replace_section  # only in the `wandb` extra

    data = read_yaml(campaign)
    ledger = read_yaml(verdicts)
    typer.echo(
        replace_section(
            url,
            RESULTS,
            results_markdown(data, ledger),
            panels=table_panels(entity, project),
            collapsed={RUN_LOG: run_log_markdown(ledger)},
        )
    )


if __name__ == "__main__":
    app()
