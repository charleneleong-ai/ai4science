"""Mirror the campaign's verdict ledger and result tables into W&B.

One run per *experiment*, not per Boltz invocation: the question a reader brings to this project
is "what was tried, what did it settle, and what is still open", and that is one row per entry in
`docs/experiments/pmhc/verdicts.yaml`. Each run carries its gate and result as notes, its verdict
as a tag, and its headline numbers in `run.summary`, so the runs table answers the question
without opening a single run.

The measurement tables ride along on the runs they belong to as `wandb.Table`s, and the host is
taken from the API rather than assumed to be wandb.ai.

    uv run --extra wandb python scripts/log_wandb.py --entity chaleong --project pmhc-g12d-binder
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import typer
from build_wandb_report import api_overrides, read_yaml

app = typer.Typer(add_completion=False)

# Which committed table belongs to which experiment, so a run carries its own evidence rather
# than every run carrying all of it.
TABLES: dict[str, tuple[str, Path]] = {
    "cohort-replicate": ("ranking", Path("data/cofold_spread_pooled.csv")),
    "per-design-denominator": ("replicates", Path("data/cofold_replicates.csv")),
    "wildtype-finalists": ("tier1", Path("data/cofold_wildtype.csv")),
}

# The numbers worth putting in the runs table itself, by experiment.
SUMMARIES: dict[str, dict[str, str]] = {
    "cohort-replicate": {
        "best_drop": "headline.best_drop",
        "bar": "headline.bar",
        "above_bar": "headline.above_bar",
        "tied_bar": "headline.tied_bar",
        "below_bar": "headline.below_bar",
        "clear_zero": "headline.clear_zero",
        "screened": "headline.screened",
        "novelty_admissible": "headline.novelty_admissible",
    },
    "bar-reproducibility": {
        "bar": "bar.pooled",
        "bar_invocations": "bar.invocations",
        "bar_clean_9": "bar.clean_9",
        "bar_spread": "bar.pooled_spread",
    },
    "per-design-denominator": {
        "between_sd": "variance.designs.between_sd",
        "within_sem": "variance.designs.within_sem",
        "run_effect": "variance.designs.run_effect",
    },
    "tier1-resolvability": {
        "control_wt_bias": "headline.tier1_control_bias",
        "control_wt_bias_sem": "headline.tier1_control_bias_sem",
        "invocations": "headline.tier1_control_invocations",
    },
    "wildtype-finalists": {
        "control_bar": "tier1_paired.control_bar",
        "control_bar_sem": "tier1_paired.control_bar_sem",
        "cohort_mean": "tier1_paired.cohort_mean",
        "invocations": "tier1_paired.invocations",
        "predictions": "tier1_paired.predictions",
        "below_bar": "tier1_paired.below_bar",
        "tied_bar": "tier1_paired.tied_bar",
        "above_bar": "tier1_paired.above_bar",
        "bulge_core_fisher_p": "tier1_paired.paired_result.bulge_core_fisher_p",
        "bulge_core_fisher_p_unpaired": "tier1_paired.withdrawn_claim.bulge_core_fisher_p",
    },
    "escape-cone": {
        "escape_cone_pct_at_8A": "escape_cone.by_reach_angstrom.8",
        "free_rays_at_8A": "escape_cone.free_rays_at_8A",
        "rays": "escape_cone.rays",
    },
}


def dig(data: dict[str, Any], path: str) -> Any:
    """`a.b.8` into nested mappings, tolerating integer-looking keys."""
    node: Any = data
    for part in path.split("."):
        node = node[int(part)] if part.isdigit() and int(part) in node else node[part]
    return node


def table_of(path: Path, wandb: Any) -> Any | None:
    if not path.exists():
        return None
    with path.open() as handle:
        rows = list(csv.reader(handle))
    return wandb.Table(columns=rows[0], data=rows[1:]) if len(rows) > 1 else None


@app.command()
def main(
    entity: str = typer.Option(..., help="W&B entity"),
    project: str = typer.Option("pmhc-g12d-binder"),
    campaign: Path = typer.Option(Path("docs/experiments/pmhc/campaign.yaml")),
    verdicts: Path = typer.Option(Path("docs/experiments/pmhc/verdicts.yaml")),
    dry_run: bool = typer.Option(False, help="print what would be logged and exit"),
) -> None:
    data = read_yaml(campaign)
    ledger = read_yaml(verdicts)

    planned = []
    for experiment in ledger["experiments"]:
        name = experiment["id"]
        summary = {
            key: dig(data, path) for key, path in SUMMARIES.get(name, {}).items()
        }
        table = TABLES.get(name)
        planned.append((experiment, summary, table if table and table[1].exists() else None))

    if dry_run:
        for experiment, summary, table in planned:
            typer.echo(
                f"{experiment['id']:26s} {experiment['verdict']:9s}"
                f" summary={list(summary) or '-'} table={table[1].name if table else 'MISSING'}"
            )
        typer.echo(f"\n{len(planned)} runs would be written to {entity}/{project}")
        return

    import wandb  # heavy, and only needed when actually logging

    host = api_overrides().get("base_url")
    for experiment, summary, table in planned:
        run = wandb.init(
            entity=entity,
            project=project,
            name=experiment["id"],
            notes=f"{experiment['what']}\n\nGate: {experiment['gate']}\n\n"
            f"{str(experiment['result']).strip()}\n\n"
            f"Consequence: {str(experiment.get('consequence', '')).strip()}",
            tags=["pmhc", "challenge5", f"verdict:{experiment['verdict']}"],
            config={"gate": experiment["gate"], "pr": experiment.get("pr")},
            reinit=True,
            settings=wandb.Settings(base_url=host) if host else None,
        )
        run.summary.update({"verdict": experiment["verdict"], **summary})
        if table and (logged := table_of(table[1], wandb)) is not None:
            run.log({f"table/{table[0]}": logged})
        typer.echo(f"  {experiment['id']:26s} {run.url}")
        run.finish()

    typer.echo(
        f"\n{len(planned)} runs in {entity}/{project}."
        " What is not established is in the project README, not in a run:\n  "
        + "\n  ".join(f"- {str(item).strip()}" for item in ledger["not_established"])
    )


if __name__ == "__main__":
    app()
