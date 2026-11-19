"""Sync the VQCC verdict ledger to W&B: one run per experiment, then the project report.

Every row of `docs/experiments/verdicts.yaml` becomes a run whose id is the experiment name, so a
re-sync updates it in place. Submissions carry their per-task costs and ratios from the matching
parts record, priced against `DEFAULT_CALIBRATION`. A `ledger` run holds the verdict table and
the committed progress chart. The report is created once (its URL kept in `wandb-report.json`);
after that Highlights, Results and the Full run log are regenerated from the ledger and the
newest leaderboard snapshot, while Task, Approach, What we learned and Next stay hand-written.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import cached_property
from typing import Any

import typer
import wandb
import wandb_workspaces.reports.v2 as wr
import yaml
from autoresearch.wandb_report import replace_section, results_table, verdict_table
from build_artifact import (
    EXPERIMENTS_DIR,
    GOLD_RANK,
    LEADERBOARD_PATH,
    PROGRESS_PATH,
    PR_URL,
    STATUS_DATE,
    VERDICTS_PATH,
    Standing,
    parts_path_for,
    scored_submissions,
)

from quantum_compilation.cost import DEFAULT_CALIBRATION, TASKS, task_ratio

REPORT_PATH = EXPERIMENTS_DIR / "wandb-report.json"
ARTIFACT_URL = "https://claude.ai/artifact/EzjSoAxvhTvi9pWgdJDEDq"
FAMILIES = {"g": "gate", "v": "submission", "r": "research"}
LEDGER_TABLE = "table/ledger"
PROGRESS_KEY = "progress"
# Hand-written: refresh at each milestone, beside the verdict that changes them.
TAKEAWAYS = (
    "- Every score so far came from picking newer public parts, not from our own compression.\n"
    "- The local verifier matched the leaderboard exactly three times; v4's -1,527 ppm is one "
    "failed `dqi` shot, the step measured in g1.\n"
    "- Our one compression attempt (r7, `invsqrt` lookup redesign) was closed on cost."
)


def first_sentence(text: str) -> str:
    plain = re.sub(r"[*`]", "", " ".join(text.split()))
    return re.split(r"(?<=[.!?])\s", plain, maxsplit=1)[0]


@dataclass
class LedgerSync:
    entity: str
    project: str

    @cached_property
    def ledger(self) -> dict[str, Any]:
        return yaml.safe_load(VERDICTS_PATH.read_text())

    @cached_property
    def submissions(self) -> list[dict[str, Any]]:
        return scored_submissions(self.ledger)

    def task_costs(self, entry: dict[str, Any]) -> dict[str, float]:
        path = parts_path_for(str(entry.get("submission", "")))
        if not path.is_file():
            return {}
        parts = yaml.safe_load(path.read_text())["parts"]
        return {task: float(parts[task]["our_measurement"]["cost"]) for task in TASKS}

    def summary(self, entry: dict[str, Any]) -> dict[str, Any]:
        out: dict[str, Any] = {"verdict": entry["verdict"]}
        if entry.get("score") is not None:
            out["score"] = float(entry["score"])
        for task, cost in self.task_costs(entry).items():
            out[f"cost/{task}"] = cost
            out[f"ratio/{task}"] = task_ratio(task, cost, DEFAULT_CALIBRATION)
        return out

    def init(self, run_id: str, **kwargs: Any) -> Any:
        return wandb.init(
            entity=self.entity, project=self.project, id=run_id, name=run_id, resume="allow",
            **kwargs,
        )

    def sync_entry(self, entry: dict[str, Any]) -> dict[str, Any]:
        family = FAMILIES.get(entry["name"][0], "other")
        run = self.init(
            entry["name"],
            job_type=family,
            group=family,
            tags=[family, f"verdict:{entry['verdict']}"],
            notes=f"{first_sentence(entry['result'])} PR #{entry['pr']}.",
            config={k: entry[k] for k in ("date", "baseline", "gate", "pr", "submission") if k in entry},
        )
        summary = self.summary(entry)
        run.summary.update(summary)
        url = run.url
        run.finish()
        return {"tag": entry["name"], "wandb_url": url, "date": entry["date"], "pr": entry["pr"],
                **summary}

    def sync_ledger(self, rows: list[dict[str, Any]]) -> None:
        run = self.init("ledger", job_type="ledger", tags=["ledger"],
                        notes=f"Verdict ledger and progress chart; dashboard {ARTIFACT_URL}")
        table = wandb.Table(
            columns=["experiment", "date", "verdict", "score", "pr"],
            data=[[r["tag"], r["date"], r["verdict"], r.get("score"), r["pr"]] for r in rows],
        )
        run.log({LEDGER_TABLE: table, PROGRESS_KEY: wandb.Image(str(PROGRESS_PATH))})
        run.finish()

    def highlights(self) -> str:
        board = json.loads(LEADERBOARD_PATH.read_text())
        latest = self.submissions[-1]
        standing = Standing(board["scores"], float(latest["score"]))
        scores = board["scores"]
        lines = [
            f"| {e['name']} | {e['date']} | {e['submission']} | "
            + (f"**{e['score']}**" if e is latest else f"{e['score']}") + " |"
            for e in self.submissions
        ]
        return (
            f"## Highlights\n\n*Status {STATUS_DATE}.* Leaderboard **{latest['score']}** "
            f"(submission {latest['submission']}), rank {standing.best:,}-{standing.worst:,} of "
            f"{standing.teams:,}; leader {standing.leader:.3f}, gold line (rank {GOLD_RANK}) "
            f"{scores[GOLD_RANK - 1]:.3f}, median {scores[standing.teams // 2]:.3f} "
            f"(`{LEADERBOARD_PATH.name}`).\n\n"
            "| Submission | Date | Ref | Leaderboard |\n|---|---|---|---|\n"
            + "\n".join(lines)
            + f"\n\n{TAKEAWAYS}\n\nDashboard: [VQCC compression log]({ARTIFACT_URL}). "
            f"PR: [#105]({PR_URL}/105)."
        )

    def results(self, rows: list[dict[str, Any]]) -> str:
        scored = [r for r in rows if "score" in r]
        gates = "\n".join(f"- **{g}**: {text}" for g, text in self.ledger["gates"].items())
        entries = [
            {**e, "label": e["name"], "result": first_sentence(e["result"])}
            for e in self.ledger["experiments"]
        ]
        return (
            "## Results\n\n### Leaderboard submissions\n\n"
            + results_table(scored, ["date", "score"] + [f"ratio/{t}" for t in TASKS],
                            label="Submission")
            + "\n### Verdicts\n\nGates:\n\n" + gates + "\n\n"
            + verdict_table(entries, PR_URL)
        )


def panels() -> list[Any]:
    return [
        wr.BarPlot(title="Leaderboard score", metrics=[wr.SummaryMetric("score")], orientation="h"),
        wr.WeavePanelSummaryTable(table_name=LEDGER_TABLE),
        wr.MediaBrowser(media_keys=[PROGRESS_KEY]),
    ]


def narrative() -> list[Any]:
    """The report's skeleton on first creation; Highlights and Results are filled in after."""
    return [
        wr.TableOfContents(),
        wr.H2(text="Highlights"),
        wr.H2(text="Task"),
        wr.MarkdownBlock(text=(
            "Kaggle `google-verifiable-quantum-circuit-compilation`. The score is the geometric mean "
            "over `invsqrt`, `dqi` and `dfthcsos` of reference cost / submitted cost, each cost "
            "measured by the task's shipped verifier."
        )),
        wr.H2(text="Approach"),
        wr.MarkdownBlock(text=(
            "A local cost model read off the three shipped verifiers (`quantum_compilation.cost`), "
            "calibrated on the reference submission (gate G0). Each candidate is SHA-256 checked "
            "against its source manifest, scored locally with all three verifiers, and submitted "
            "only when it clears the gates below. Provenance per submission is in "
            "`docs/experiments/<submission>-parts.yaml`."
        )),
        wr.H2(text="Results"),
        wr.H2(text="What we learned"),
        wr.MarkdownBlock(text=(
            "| Lever | Measured gain |\n|---|---|\n"
            "| Newer public parts (v2 -> v3 -> v4) | x23.18, x1.094, x1.111 |\n\n"
            "**Closed / do not retry**\n\n"
            "- THC factorization refit (r1): lambda_eff floor 46.786 leaves 3.1% of one factor.\n"
            "- brot 7 or 9, shape re-splits (r2): accuracy cap breached, or width doubles.\n"
            "- Cheaper J=5 loader synthesis (r4): the public planner is already at generic rank.\n"
            "- `invsqrt` lookup redesign (r7): synthesis 515-517 Toffolis against a 231 budget."
        )),
        wr.H2(text="Next"),
        wr.MarkdownBlock(text=(
            "- `dqi`: the small open patches (closing Toffolis, zero-word guard, CZ commuting).\n"
            "- Milestone snapshot 2026-11-20 23:59 UTC; final deadline 2027-01-20 23:59 UTC; 5 "
            "submissions a day."
        )),
    ]


def create_report(sync: LedgerSync) -> str:
    report = wr.Report(
        entity=sync.entity, project=sync.project, title="VQCC experiment log",
        description="Verifiable quantum circuit compilation: submissions, verdicts and what is closed.",
    )
    report.blocks = narrative()
    report.save()
    return report.url


def main(
    entity: str = typer.Option(None, help="W&B entity; defaults to the logged-in one"),
    project: str = typer.Option("vqcc"),
) -> None:
    sync = LedgerSync(entity or wandb.Api().default_entity, project)
    rows = [sync.sync_entry(e) for e in sync.ledger["experiments"]]
    sync.sync_ledger(rows)
    url = json.loads(REPORT_PATH.read_text())["url"] if REPORT_PATH.exists() else create_report(sync)
    replace_section(url, "Highlights", sync.highlights())
    url = replace_section(
        url, "Results", sync.results(rows), panels=panels(),
        collapsed={"Full run log": results_table(rows, ["date", "verdict", "score", "pr"], label="Run")},
    )
    REPORT_PATH.write_text(json.dumps({"url": url}, indent=2) + "\n")
    typer.echo(f"{len(rows)} runs synced; report {url}")


if __name__ == "__main__":
    typer.run(main)
