#!/usr/bin/env -S uv run --quiet --with wandb-workspaces --with typer --with rich python
"""Build the campaign's W&B report with the run's own panels embedded, not just prose.

A markdown-only report restates results the reader then has to go and find. This one pulls
the logged tables and the folded complexes into panel grids bound to the run, so every claim
sits next to the measurement it came from.

Panel keys are read back off the run rather than hardcoded, so a key renamed in tracking.py
surfaces here as a missing panel instead of a silently empty one.
"""

from __future__ import annotations

import typer
import wandb
import wandb_workspaces.reports.v2 as wr
from rich.console import Console

app = typer.Typer(add_completion=False, help="Publish the campaign report to Weights & Biases.")

NARRATIVE_INTRO = """\
Design a TNF-α binder that **binds at pH 7.4 and releases at pH 6.0**, **cross-reacts with
mouse TNF-α**, and binds tightly — ranked in that order. The organisers' FAQ is explicit that a
weak but clearly pH-selective binder beats a high-affinity non-selective one, so selectivity is
optimised first and affinity last.

The epitope was computed, not cited: TNF-α residues within 4.5 Å of TNFR2 in **3ALQ**,
intersected with a human/mouse alignment (**124/157 = 79.0% identity**), give **33 contact
positions — 24 conserved, 9 divergent**. The nine divergent positions are the cross-reactivity
trap. **R31→Q** carries the largest solvent-exposed area in the epitope (193 Å²) and sits beside
R32 in the 29–33 loop that both adalimumab and infliximab engage, which is precisely why neither
is mouse cross-reactive. Those nine became coldspots; two fully-conserved sub-patches were
targeted instead."""

NARRATIVE_FUNNEL = """\
## The funnel is the result

**57 folded designs → 3 submitted → 2 passed Proteinbase novelty**, a yield of 3.5%. Reporting
the 3 alone would be marketing; the denominators are the finding. Composition is where the
generator actually failed — 46 of the 53 designs surviving a memorisation check were
low-complexity alanine or glycine runs that fold confidently and mean nothing."""

NARRATIVE_LEVERS = """\
## Only one design lever moved the interface

Three of the four levers I expected to matter did nothing. Modelling the target as its true
homotrimer was the only win, and truncating it to a dimer silently cost 0.11."""

NARRATIVE_CONFIDENCE = """\
## Global ipTM is unusable on a multi-chain target

TNF-α is a homotrimer, so three of the six chain pairs in a binder/trimer complex are
target-internal and score 0.85–0.96 regardless of the binder, dominating the global average.
**Measured across all 57 designs, global ipTM overstates the designed interface by +0.405** —
roughly half its own value. Every number here is instead a `pair_chains_iptm` entry between the
binder chain and each target chain, symmetrised over both matrix directions (which differ by up
to 0.08), averaged across target chains, then averaged across diffusion samples.

Aggregation is by **mean, not best-of-n**. Ranking on the best sample selects for sampling
variance: best and mean differ by up to 0.35 within a single design, and the ordering changes.

### Confidence was not sufficiency

Six designs cleared the ~0.50 bar for a credible de novo interface and **none was submittable**.
Three were 59% alanine or 37% glycine. The clean ones — `pocket_free_1` and `pocket_free_2` —
had real sequences and real confidence, and put only **13–18% of their interface on the epitope
they were conditioned on**. They were binding TNF-α somewhere else. Epitope focus, not
confidence, was the operative constraint, and a gate reading confidence alone passes all six."""

NARRATIVE_PH = """\
## The pH switch, and three refuted mechanisms

PROPKA 3.5.1 over PyRosetta-built side chains, scanning **200 histidine substitutions**: 18 land
in the ideal 6.3–7.2 pKa window with the histidine at the interface (ΔSASA ≥ 10 Å²), 17 more in
the usable 6.0–7.6 band. Multi-histidine combinations shift each other's pKa by **≤ 0.21**, so
cooperative switches are viable rather than self-defeating.

Three mechanisms I reasoned my way to were all wrong, and the scan caught each: **burying the
histidine** lowers its pKa rather than raising it (r = −0.37; buried sites land at 3.1–3.6);
**a nearby cation** has the wrong sign (r = +0.11); **carboxylate proximity** gives r = +0.21
with counter-examples at net +6 charge sitting at pKa 3.6. A brute-force scan plus a physics
solver beat every geometric rule I invented."""

NARRATIVE_NOVELTY = """\
## Fold, not sequence, is the novelty constraint

Proteinbase scores novelty 1–4 and requires ≥3. Of the three submitted designs, the one that
**failed** had *lower* sequence identity (37.5%) than the two that passed (38.1%). What sank it
was fold: it matched **CATH 2.60.40**, the immunoglobulin-like β-sandwich, at TM 0.893 over 88.8%
coverage. The two α-bundles passed. The rule that follows — avoid Ig-like β-sandwiches, α-bundles
reach score 3 — also removed `mixed_topology` from the later BindCraft2 campaigns.

One calibration claim had to be **retracted**. I reported that our fold-novelty gate reproduced
Proteinbase's verdict 3/3; it only scored the two passes correctly because a missing identity
value defaulted to 0. With the real 38.1% the gate returns 2 where Proteinbase passed. The gate
is **conservative, not calibrated** — decisive on failure, advisory on borderline — and a
regression test now pins that."""

NARRATIVE_CLOSE = """\
## What is not established

**There is no measured binding of any kind.** No expression, no SEC, no SPR or BLI. pH
selectivity rests on predicted pKa, and no folding model represents protonation state. Mouse
cross-reactivity rests on a predicted mouse structure and conserved-epitope coverage, not on
measured binding. Interface confidence is model-dependent and unresolved, so a range is reported
rather than a number: AlphaFold2 rated its own optimised designs 0.92 where Boltz-2 rated the
same structures 0.21.

A generated sequence that has passed computational checks is a candidate. A result is a measured
biological endpoint on a preserved funnel, with controls. This produced the former."""


def media_keys(run: "wandb.apis.public.Run", prefix: str = "structure/") -> list[str]:
    """Logged media keys, read back off the run so a rename shows up as a missing panel."""
    return sorted(k for k in run.summary.keys() if k.startswith(prefix))


def table_keys(run: "wandb.apis.public.Run", prefix: str = "table/") -> list[str]:
    return sorted(k for k in run.summary.keys() if k.startswith(prefix))


@app.command()
def main(
    run_path: str = typer.Argument(..., help="entity/project/run_id of the campaign run"),
    title: str = typer.Option("Conditional TNF-α Binder — Adaptyv × Anthropic Challenge 02",
                              "--title"),
) -> None:
    """Publish the report, embedding the named run's tables and structures as panels."""
    console = Console()
    entity, project, run_id = run_path.split("/")
    run = wandb.Api().run(run_path)
    structures, tables = media_keys(run), table_keys(run)
    console.print(f"run [cyan]{run.name}[/] · {len(tables)} tables · "
                  f"{len(structures)} structures")

    runset = wr.Runset(entity=entity, project=project, name="campaign",
                       filters=f"name == '{run_id}'")

    def grid(*panels: object) -> wr.PanelGrid:
        return wr.PanelGrid(runsets=[runset], panels=list(panels))

    def summary_table(key: str) -> wr.WeavePanelSummaryTable:
        return wr.WeavePanelSummaryTable(table_name=key)

    report = wr.Report(
        entity=entity, project=project, title=title, width="fluid",
        description="Funnel, measured design levers, and why high interface confidence did "
                    "not mean a usable design.",
    )
    report.blocks = [
        wr.MarkdownBlock(text=NARRATIVE_INTRO),
        wr.CalloutBlock(text=[
            "Every figure below is live from the campaign run — ",
            wr.Link(text=f"{project}/{run.name}", url=run.url),
            ". Structures are the folded complexes themselves; rotate them in place.",
        ]),
        wr.TableOfContents(),

        wr.MarkdownBlock(text=NARRATIVE_FUNNEL),
        grid(summary_table("table/funnel"),
             wr.BarPlot(title="Survivors per stage", metrics=["funnel/surviving"])),

        wr.MarkdownBlock(text=NARRATIVE_LEVERS),
        grid(summary_table("table/design_levers")),

        wr.MarkdownBlock(text=NARRATIVE_CONFIDENCE),
        grid(summary_table("table/designs")),
        wr.MarkdownBlock(
            text="The `inflation` column is `global_iptm − binder_to_target`: the part of the "
                 "global score contributed by the target folding against itself."),
        grid(wr.MediaBrowser(title="Folded complexes, best interface first",
                             media_keys=structures, num_columns=3)),

        wr.MarkdownBlock(text=NARRATIVE_PH),
        grid(summary_table("table/pka_scan2")),
        grid(summary_table("table/his_combo_pocket_free_0")),

        wr.MarkdownBlock(text=NARRATIVE_NOVELTY),
        wr.MarkdownBlock(text=NARRATIVE_CLOSE),
    ]
    report.save()
    console.print(f"[green]report[/] {report.url}")


if __name__ == "__main__":
    app()
