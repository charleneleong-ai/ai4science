"""W&B runs and report for the Kaggle submissions, rendered from `docs/experiments/casmi/verdicts.yaml`.

One run per scored submission (id `v<version>`, so a re-sync updates rather than duplicates) and a
`trajectory` run holding the score-by-version table. Everything is written to run summaries, never
to history: a resumed run silently drops history steps it has already passed. The report's
hand-written sections come from `report.md` beside the yaml; Results and Full run log are generated.
Re-run after each score lands.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import wandb
import wandb_workspaces.reports.v2 as wr
import yaml
from autoresearch.wandb_report import results_table, splice_section, verdict_table

PROJECT = "casmi26"
REPO = "https://github.com/charleneleong-ai/ai4science"
SAMPLE_MOLECULES = 12


def load(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text())


def annotate(subs: list[dict[str, Any]], noise: float) -> list[dict[str, Any]]:
    """Each submission with the best score before it and whether it cleared gate G1."""
    best = None
    out = []
    for s in subs:
        prior = best
        cleared = prior is None or s["score"] > prior + noise
        best = s["score"] if prior is None else max(prior, s["score"])
        out.append({**s, "best_before": prior, "gate_g1": cleared, "best_so_far": best})
    return out


def top1_table(outputs: Path, version: int) -> wandb.Table | None:
    """Top-1 structures for the first molecules of a kept run output, as 3D molecules."""
    sub = outputs / f"v{version}" / "submission.csv"
    if not sub.exists():
        return None
    table = wandb.Table(columns=["molecule_id", "top1_smiles", "structure"])
    with sub.open() as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            if i == SAMPLE_MOLECULES:
                break
            smi = row["smiles"].split(";")[0]
            table.add_data(row["molecule_id"], smi, wandb.Molecule.from_smiles(smi))
    return table


def sync_run(s: dict[str, Any], pr: int, entity: str, outputs: Path) -> str:
    notes = (f"v{s['version']}: {s['change']}. Public MRR@25 {s['score']:.3f} "
             f"(best before: {s['best_before']}). Verdict {s['verdict']}. {REPO}/pull/{pr}")
    run = wandb.init(project=PROJECT, entity=entity, id=f"v{s['version']}", name=f"v{s['version']}",
                     resume="allow", job_type="submission", notes=notes,
                     tags=[s["family"], f"verdict:{s['verdict']}"],
                     config={k: s[k] for k in ("version", "ref", "date", "commit", "change", "family")})
    run.summary.update({"public_mrr25": s["score"], "best_so_far": s["best_so_far"],
                        "gate_g1": s["gate_g1"], "verdict": s["verdict"],
                        "commit_url": f"{REPO}/commit/{s['commit']}"})
    if (table := top1_table(outputs, s["version"])) is not None:
        run.summary["top1_structures"] = table
    url = run.url
    run.finish()
    return url


def sync_trajectory(subs: list[dict[str, Any]], entity: str) -> None:
    run = wandb.init(project=PROJECT, entity=entity, id="trajectory", name="trajectory",
                     resume="allow", job_type="summary", tags=["trajectory"],
                     notes="Public MRR@25 and best-so-far by notebook version, one row per scored submission.")
    run.summary["trajectory"] = wandb.Table(
        columns=["version", "public_mrr25", "best_so_far", "verdict"],
        data=[[s["version"], s["score"], s["best_so_far"], s["verdict"]] for s in subs])
    run.finish()


def narrative(path: Path) -> dict[str, str]:
    """`## Heading` -> body, from the hand-written report.md."""
    parts = re.split(r"^## (.+)$", path.read_text(), flags=re.M)[1:]
    return {h.strip(): f"## {h.strip()}\n{body.strip()}\n" for h, body in zip(parts[::2], parts[1::2])}


def results_markdown(subs: list[dict[str, Any]], pr: int, noise: float) -> str:
    rows = [{"tag": f"v{s['version']}", "wandb_url": s.get("url"), "score": s["score"],
             "verdict": s["verdict"]} for s in subs]
    entries = [{"label": f"v{s['version']}", "baseline": "—" if s["best_before"] is None else f"{s['best_before']:.3f}",
                "gate": "G1", "result": f"{s['score']:.3f}: {s['change']}", "verdict": s["verdict"], "pr": pr}
               for s in subs]
    return ("## Results\n"
            f"Public leaderboard MRR@25 for every scored submission. G1 = beats the best so far by more than "
            f"the run-to-run noise (±{noise}).\n\n"
            + results_table(rows, ["score", "verdict"], label="Version", fmt="{:.3f}")
            + "\n### Experiment verdicts\n" + verdict_table(entries, pr_url=f"{REPO}/pull"))


def run_log(subs: list[dict[str, Any]]) -> str:
    return "\n".join(
        f"- [v{s['version']}]({s.get('url', '')}) · {s['date']} · [`{s['commit']}`]({REPO}/commit/{s['commit']}) · "
        f"submission {s['ref']} · {s['score']:.3f} · {s['verdict']} — {s['change']}" for s in subs)


def section_blocks(heading: str, md: str, subs: list[dict[str, Any]], data: dict[str, Any],
                   entity: str) -> list[Any]:
    """Report blocks for one section: the hand-written markdown, or the generated content."""
    if heading == "Results":
        grid = wr.PanelGrid(runsets=[wr.Runset(entity=entity, project=PROJECT)], panels=[
            wr.BarPlot(title="Public MRR@25 by submission", metrics=["public_mrr25"]),
            wr.WeavePanelSummaryTable(table_name="trajectory")])
        return [wr.MarkdownBlock(text=results_markdown(subs, data["pr"], data["noise"])), grid]
    if heading == "Full run log":
        return [wr.MarkdownBlock(text="## Full run log\nEvery submission, oldest first.\n"),
                wr.H2(text="Runs", collapsed_blocks=[wr.MarkdownBlock(text=run_log(subs))])]
    return [wr.MarkdownBlock(text=md)]


def open_report(verdicts_path: Path, data: dict[str, Any], entity: str, sections: dict[str, str]) -> Any:
    """The project report; created from the narrative on the first sync and its URL recorded."""
    if url := data.get("report_url"):
        return wr.Report.from_url(url)
    report = wr.Report(project=PROJECT, entity=entity,
                       title="CASMI 2026 — molecule identification from tandem mass spectra",
                       description="One change per submission, scored only on the public leaderboard.")
    report.blocks = [wr.TableOfContents()] + [wr.MarkdownBlock(text=md) for md in sections.values()]
    report.save()
    verdicts_path.write_text(verdicts_path.read_text().rstrip() + f"\nreport_url: {report.url}\n")
    return report


def sync(verdicts_path: Path, outputs: Path) -> str:
    data = load(verdicts_path)
    entity = wandb.Api().default_entity
    subs = annotate(data["submissions"], data["noise"])
    for s in subs:
        s["url"] = sync_run(s, data["pr"], entity, outputs)
    sync_trajectory(subs, entity)
    sections = narrative(verdicts_path.parent / "report.md")
    report = open_report(verdicts_path, data, entity, sections)
    blocks = report.blocks
    for heading, md in sections.items():
        blocks = splice_section(blocks, heading, section_blocks(heading, md, subs, data, entity), wr,
                                drop=("Runs",) if heading == "Full run log" else ())
    report.blocks = blocks
    report.save()
    return report.url
