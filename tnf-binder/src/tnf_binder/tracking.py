"""Logging the campaign to Weights & Biases: the funnel, the confidence tables, the structures.

What is logged is chosen to make the campaign auditable by someone who was not in the session:
the elimination funnel with denominators, every design's confidence (mean across diffusion
samples, with the global-ipTM inflation shown beside it), the pKa scan the pH switch was
chosen from, and the folded complexes themselves as 3D media.

`wandb` is imported lazily -- it is an optional extra, and the rest of the package must stay
importable without it.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any

from tnf_binder.confidence import DesignConfidence, mean_or_zero, rank_designs

# (stage, candidates entering, surviving) -- the denominators are the result, not the 3 winners.
FUNNEL: tuple[tuple[str, int, int], ...] = (
    ("generated", 57, 57),
    ("not memorised (ubiquitin / VHH framework)", 57, 53),
    ("composition (Ala 59-76% rejected)", 53, 7),
    ("interface (clash / tiny / off-epitope)", 7, 3),
    ("submitted to Proteinbase", 3, 3),
    ("passed novelty >= 3/4", 3, 2),
)

CREDIBLE_INTERFACE = 0.5  # the bar a de novo binder interface must clear

# The second submission's own funnel: a His series over the chassis that survived the campaign's
# gates. Separate from FUNNEL because the denominators are different objects -- parents and
# variants, not the original 57 constructs.
SUBMISSION2_FUNNEL: tuple[tuple[str, int, int], ...] = (
    ("folded parent designs", 16, 16),
    ("fold similarity only moderate", 16, 9),
    ("composition pass", 9, 2),
    ("His variants built and folded", 15, 15),
    ("cross-reactivity >= 0.6", 15, 12),
    ("submitted", 12, 12),
)

# Measured effect of each design lever on binder-to-target ipTM. Only one of the four moved.
LEVERS: tuple[tuple[str, str, float], ...] = (
    ("target stoichiometry", "dimer -> trimer", +0.11),
    ("architecture", "de novo -> nanobody", 0.00),
    ("topology", "free -> helix-biased", -0.03),
    ("epitope", "pocket -> lobe", 0.00),
)

CAMPAIGN_CONFIG: dict[str, Any] = {
    "target": "human TNF-alpha homotrimer (P01375), TNFR2 epitope from 3ALQ",
    "human_mouse_identity": 0.790,
    "epitope_contacts": 33,
    "epitope_conserved": 24,
    "objectives": ["pH-selective release (6.0 vs 7.4)", "mouse cross-reactivity", "affinity"],
    "generators": ["BoltzGen 0.1.4", "BindCraft 1.5.3", "BindCraft2"],
    "folding_oracles": ["Boltz-2 2.2.1", "AlphaFold2 (BindCraft)"],
    "ranking_statistic": "mean binder-to-target pair ipTM across diffusion samples",
}



def funnel_surviving(stage: str, funnel: tuple[tuple[str, int, int], ...] = FUNNEL) -> int:
    """Survivors at the first funnel stage whose name contains `stage`.

    Keeps the run summary derived from FUNNEL rather than restating its numbers, which could
    otherwise drift apart silently.
    """
    for name, _entering, surviving in funnel:
        if stage in name:
            return surviving
    raise KeyError(f"no funnel stage matching {stage!r}")


@dataclass(frozen=True)
class CampaignLogger:
    """One W&B run holding the whole TNF-alpha campaign."""

    assets: Path
    predictions_dir: str = "boltz_camp"
    funnel: tuple[tuple[str, int, int], ...] = FUNNEL
    project: str = "tnf-conditional-binder"
    entity: str | None = None
    max_molecules: int = 12

    @cached_property
    def designs(self) -> list[DesignConfidence]:
        """Every folded design, best-ranked first by mean binder-to-target ipTM."""
        return rank_designs(self.assets / self.predictions_dir)

    @cached_property
    def structures(self) -> list[Path]:
        """Folded complexes worth attaching as 3D media, ranked designs first."""
        paths: list[Path] = []
        for design in self.designs:
            paths += sorted(
                (self.assets / self.predictions_dir / design.name).glob("*_model_0.cif"))
        paths += sorted(self.assets.rglob("*.pdb"))
        return paths[:self.max_molecules]

    @staticmethod
    def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
        with path.open() as handle:
            rows = list(csv.reader(handle))
        return (rows[0], rows[1:]) if rows else ([], [])

    @staticmethod
    def molecule_key(path: Path) -> str:
        """One CIF per design directory, so the directory names it; loose PDBs name themselves."""
        return path.parent.name if path.suffix == ".cif" else path.stem

    def summary(self) -> dict[str, float]:
        designs = self.designs
        summary = {
            "designs_folded": len(designs),
            "designs_submitted": funnel_surviving("submitted", self.funnel),
            "best_binder_to_target_iptm": round(designs[0].binder_to_target, 3) if designs else 0,
            "mean_global_iptm_inflation": round(mean_or_zero(d.inflation for d in designs), 3),
            "credible_interface_bar": CREDIBLE_INTERFACE,
        }
        # only the campaign funnel reaches a novelty verdict; a submission still in flight has
        # no such number, and reporting one from the wrong funnel would be worse than omitting it
        if any("novelty" in stage for stage, _, _ in self.funnel):
            summary["designs_passing_novelty"] = funnel_surviving("novelty", self.funnel)
        return summary

    def log(self, *, run_name: str | None = None, notes: str | None = None) -> str:
        """Create the run, log everything, return its URL."""
        import wandb

        designs = self.designs
        run = wandb.init(
            project=self.project, entity=self.entity, name=run_name, notes=notes,
            job_type="campaign", tags=["protein-design", "tnf-alpha", "adaptyv", "ph-switch"],
            config=CAMPAIGN_CONFIG,
        )

        for stage, entering, surviving in self.funnel:
            wandb.log({"funnel/entering": entering, "funnel/surviving": surviving,
                       "funnel/stage": stage})

        tables = {
            "funnel": wandb.Table(columns=["stage", "entering", "surviving", "attrition"],
                                  data=[[s, e, v, round(1 - v / e, 3)]
                                        for s, e, v in self.funnel]),
            "design_levers": wandb.Table(columns=["lever", "change", "delta_iptm"],
                                         data=[list(row) for row in LEVERS]),
            "designs": wandb.Table(
                columns=["design", "samples", "binder_to_target_mean", "binder_to_target_best",
                         "global_iptm", "target_internal", "inflation", "binder_ptm"],
                data=[[d.name, len(d.samples), round(d.binder_to_target, 3),
                       round(d.binder_to_target_best, 3), round(d.global_iptm, 3),
                       round(d.target_internal, 3), round(d.inflation, 3),
                       round(d.binder_ptm, 3)] for d in designs]),
        }
        for path in sorted(self.assets.glob("*.csv")):
            columns, rows = self.read_csv(path)
            if columns:
                tables[path.stem] = wandb.Table(columns=columns, data=rows)

        molecules = {f"structure/{self.molecule_key(p)}": wandb.Molecule(str(p))
                     for p in self.structures}

        wandb.log({f"table/{k}": v for k, v in tables.items()} | molecules)
        wandb.summary.update(self.summary())
        url = run.url
        run.finish()
        return url
