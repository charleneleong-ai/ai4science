"""Choose the designs to submit, on the objectives in the order the brief ranks them.

The brief ranks pH-selectivity first, mouse cross-reactivity second, affinity third. Two further
constraints are not objectives but gates, and both can end a slot before anything is measured:

- **Novelty.** Proteinbase scores every submission 1-4 and participants report rejections at 2.
  "De novo" (4) needs sequence identity <=30% *and* less than moderate structural similarity; 3
  needs one or the other. This is why the 3DS18 chassis is not submitted despite carrying the
  strongest switch measured here: MMseqs2 finds its own crystal, 6wrv_C, at 100% identity over
  full coverage, and a three-histidine variant of it at 96.7%. That is a 1, whatever it does at
  pH 6.0. The scale is applied here as a conservative local proxy -- the fold half measures more
  similarity than Proteinbase's own pipeline does, so a `high` verdict reliably predicts
  rejection while a `moderate` one may still pass -- and a sequence hit is only counted as
  homology when it is statistically significant, because 58% identity over 17 residues at an
  E-value of 48 is not a relative.
- **Diversity.** Selection for the shared plate is made from roughly 375 designs per track by a
  workflow judging design quality, novelty and method novelty, and it will not take twenty
  near-identical sequences. Slots are therefore spread across chassis, and two designs from one
  chassis must differ by at least `MIN_SITE_DIFFERENCE` substitutions.

How the switch is scored, and why it is not the verifier's own ranking column. The de novo
chassis has one usable structure per design, so `ph_switch.verdict` returns
`single structure - not replicated` for every position and every field derived from it
collapses: `both_releasing` is False on all 640 rows of the agreement table, `n_releasing` is 0,
`total_release_consistent` is 0 and `fold_weaker_at_ph6` is 1.0 for every design -- while 359 of
those rows have a linkage past the coupling floor on both methods. Read as a ranking, the
verifier's own columns would report a false negative for all 256 variants.

The replication axis that *is* available on one pose is the two independent pKa methods, so a
position counts here only when PROPKA and pKAI both put its linkage past the floor in the
releasing direction, and it contributes the **smaller** of the two methods' log-unit estimates.
A design's switch is the sum over those positions. Nothing is ranked on the better of the two
methods, and nothing is ranked on repacks of a single pose.

Two definitions of "the switch" therefore coexist and disagree: `ph_switch` prints "0 of N
positions are coupled" for the same designs this table scores non-zero. TODO: the deeper fix is
to give the verifier a cross-method replication axis -- `ph_switch.agreement` already groups on
exactly that -- so the statistic lives with the measurement, and to make `per_design` refuse to
emit an unflagged 0 / 1.0 for a chassis its own rank key cannot measure. It needs the pKa pass
re-run, so it is recorded rather than done.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from csvio import group_by, number, read_rows, same_design, write_rows

# The two pKa methods `ph_switch` writes a delta_/release_ column pair for. Named rather than
# discovered by key prefix, so that an added or renamed column cannot silently change the score
# and so that "both methods agree" is a stated requirement rather than a length check.
METHODS = ("propka", "pkai")

# An alignment this weak is a coincidence, not a relative, however high its percent identity.
SIGNIFICANT_EVALUE = 1e-3

# Proteinbase's published bands: high similarity is >=70% coverage at TM>=0.8, moderate is >70%
# coverage at TM>=0.5, and "de novo" needs identity <=30% as well as less than moderate
# similarity. From the competition Terms, as recorded with its source at
# `glp1r/docs/strategy.md:179`; `tnf_binder.submission` records the same 3-of-4 floor and the
# 20-row cap from the live submit page, which disagreed with the FAQ and won.
HIGH_TM, MODERATE_TM, SIMILARITY_COVERAGE = 0.8, 0.5, 0.70
DE_NOVO_IDENTITY, IDENTITY_KNOWN = 30.0, 70.0
BAND_ORDER = ("less than moderate", "moderate", "high")

# Two submitted sequences from one chassis must differ by at least this many substitutions.
MIN_SITE_DIFFERENCE = 2

OUT_COLUMNS = [
    "rank", "name", "sequence", "length", "parent", "epitope", "substitutions",
    "switch_log_units", "fold_weaker_at_ph6", "coupled_sites", "epitope_conservation",
    "parent_iptm", "parent_pae", "unpaired_cys", "novelty_score", "novelty_note",
]
REJECTED_COLUMNS = ["name", "parent", "switch_log_units", "coupled_sites", "chosen", "reason"]

app = typer.Typer(add_completion=False)


def read_m8(path: Path, columns: Sequence[str]) -> list[dict[str, str]]:
    """A tabular BLAST-format hit table. The first parser of one in this repo."""
    with path.open() as handle:
        return [
            dict(zip(columns, line.rstrip("\n").split("\t"), strict=False))
            for line in handle if line.strip()
        ]


def coupling_floor(sidecar: Path | None, override: float | None) -> float:
    """The linkage floor, taken from the verifier's own resolved thresholds.

    `ph_switch` makes this a CLI option and writes the resolved value beside its outputs as
    `<out_summary>.thresholds.json`, because an earlier revision of this campaign hardcoded its
    thresholds and shipped two epitopes its own gate had rejected. Restating the number here
    would defeat that mechanism: sweeping the floor upstream would silently desync the
    selection from the measurement it is made of.
    """
    if override is not None:
        return override
    if sidecar is None or not sidecar.exists():
        raise typer.BadParameter(
            "pass --thresholds pointing at the verifier's <out_summary>.thresholds.json, or "
            "--coupling-floor to state the floor the selection used")
    return float(json.loads(sidecar.read_text())["coupling_floor"])


@dataclass(frozen=True)
class FoldHit:
    """One structural hit, in the band Proteinbase's scale would read it as."""

    target: str
    tm: float
    coverage: float

    @property
    def band(self) -> str:
        if self.coverage < SIMILARITY_COVERAGE or self.tm < MODERATE_TM:
            return "less than moderate"
        return "high" if self.tm >= HIGH_TM else "moderate"

    @property
    def rank(self) -> tuple[int, float]:
        """Band first, so coverage is never traded away for TM-score.

        Ranking on TM alone picks the wrong hit: 0.95 over 40% of the design is `less than
        moderate` on this scale and costs nothing, while 0.82 over 90% is `high` and is what a
        rejection is made of. The band decides; TM only breaks ties inside it.
        """
        return (BAND_ORDER.index(self.band), self.tm)


NO_FOLD_HIT = FoldHit("-", 0.0, 0.0)


@dataclass(frozen=True)
class Novelty:
    """One design's predicted Proteinbase band, from the sequence and the fold half."""

    identity: float
    identity_target: str
    fold: FoldHit

    @property
    def score(self) -> int:
        above_de_novo = self.identity > DE_NOVO_IDENTITY
        if self.identity > IDENTITY_KNOWN and self.fold.band in ("high", "moderate"):
            return 1
        if self.fold.band == "high" or (above_de_novo and self.fold.band == "moderate"):
            return 2
        if above_de_novo or self.fold.band == "moderate":
            return 3
        return 4

    @property
    def note(self) -> str:
        identity = (f"{self.identity:.0f}% to {self.identity_target}"
                    if self.identity else "no significant sequence hit")
        return (f"{identity}; fold {self.fold.band} "
                f"(TM {self.fold.tm:.2f} over {self.fold.coverage:.0%} to {self.fold.target})")


def novelty_of(design: str, sequence_hits: Iterable[dict[str, str]],
               fold_hits: Iterable[dict[str, str]]) -> Novelty:
    """The strongest significant sequence hit and the strongest fold band -- what a gate sees."""
    mine = [h for h in sequence_hits if same_design(design, h["query"])
            and number(h["evalue"], 1e9) <= SIGNIFICANT_EVALUE]
    best = max(mine, key=lambda h: number(h["pident"]), default=None)

    folds = [FoldHit(h["target"], number(h["qtmscore"]), number(h["qcov"]))
             for h in fold_hits if same_design(design, h["query"])]
    return Novelty(
        identity=number(best["pident"]) if best else 0.0,
        identity_target=best["target"] if best else "-",
        fold=max(folds, key=lambda hit: hit.rank, default=NO_FOLD_HIT),
    )


@dataclass(frozen=True)
class Candidate:
    """One variant, scored on the switch both methods agree about."""

    name: str
    sequence: str
    substitutions: str
    switch: float
    coupled: int

    @property
    def parent(self) -> str:
        return self.name.split("-")[0]

    @property
    def sites(self) -> frozenset[str]:
        """The substituted positions, read off the name.

        TODO: `substitutions` is already a column and `his_variants.Campaign` knows the parent
        and the site labels at write time, so both this and `parent` are reading a name used as
        a serialisation format for fields that exist -- the fix is to add them to
        `his_variants.DESIGN_COLUMNS`. Not done now for two reasons: regenerating that table
        means rebuilding the PyRosetta poses, and reading the labels from `substitutions`
        instead would make the `-wt` control parse to the empty set, putting it one
        substitution from every single-site sibling and newly suppressing it below. That moves
        a selection being made now, so it is a deliberate decision rather than a refactor.
        """
        return frozenset(self.name.split("-")[1:])

    def too_similar_to(self, other: Candidate) -> bool:
        """Whether these two are one design twice: same chassis, one substitution apart.

        The generator enumerates every subset of a design's five chosen positions, so the top of
        a switch ranking is mostly one site set and its own near-subsets, which are not
        independent designs -- the stronger one was measured, and the other is that design minus
        a site whose contribution the same table already reports. Spending two of twenty wells
        on them buys one hypothesis twice.

        The test is the symmetric difference, not containment: containment would let the maximal
        five-site set exclude every other variant of its chassis, which costs the strongest
        switch measured here all but one of its slots. Two substitutions apart is the least that
        makes two sequences a different question.
        """
        return (self.parent == other.parent
                and len(self.sites ^ other.sites) < MIN_SITE_DIFFERENCE)

    def row(self, parent: Mapping[str, str], novelty: Novelty,
            conservation: Mapping[str, float]) -> dict[str, object]:
        epitope = parent.get("epitope", "unknown")
        return {
            "name": self.name,
            "sequence": self.sequence,
            "length": len(self.sequence),
            "parent": self.parent,
            "epitope": epitope,
            "substitutions": self.substitutions,
            "switch_log_units": round(self.switch, 2),
            # The same quantity, under the same name, that `ph_switch.per_design` computes.
            "fold_weaker_at_ph6": round(10.0 ** self.switch, 1),
            "coupled_sites": self.coupled,
            "epitope_conservation": conservation.get(epitope, 0.0),
            "parent_iptm": parent.get("iptm", ""),
            "parent_pae": parent.get("pae", ""),
            "unpaired_cys": parent.get("unpaired_cys", "") or "0",
            "novelty_score": novelty.score,
            "novelty_note": novelty.note,
        }


def score_switch(rows: Iterable[dict[str, str]], floor: float) -> tuple[float, int]:
    """Sum the smaller method estimate over positions both methods call releasing."""
    total, coupled = 0.0, 0
    for row in rows:
        if not all(number(row.get(f"delta_{method}")) <= -floor for method in METHODS):
            continue
        total += min(number(row.get(f"release_{method}")) for method in METHODS)
        coupled += 1
    return total, coupled


def candidates(switch: Sequence[dict[str, str]], variants: Sequence[dict[str, str]],
               floor: float) -> list[Candidate]:
    by_design = group_by(switch, lambda row: row["design"])
    out = []
    for variant in variants:
        name = variant["name"]
        total, coupled = score_switch(by_design.get(name, []), floor)
        out.append(Candidate(name, variant["sequence"], variant.get("substitutions", ""),
                             total, coupled))
    return out


def admissible(candidate: Candidate, novelty: Mapping[str, Novelty],
               passes: Mapping[str, bool], floor: int, allow_unpaired_cys: bool,
               unpaired: Mapping[str, float]) -> str:
    """Why this candidate cannot take a slot, or an empty string if it can."""
    parent = candidate.parent
    if parent not in novelty:
        return "no novelty measurement for its parent"
    if novelty[parent].score < floor:
        return f"novelty {novelty[parent].score} below {floor} ({novelty[parent].note})"
    if not passes.get(parent, False):
        return "parent fails composition triage"
    if not allow_unpaired_cys and unpaired.get(parent, 0.0) > 0:
        return f"parent carries {unpaired[parent]:.0f} unpaired cysteine(s)"
    if candidate.coupled < 1:
        return "no position where both methods agree on release"
    return ""


def choose(ranked: Sequence[Candidate], slots: int) -> list[Candidate]:
    """Each parent's best first, then each parent's second, and so on until the slots are full.

    Filling greedily by switch instead would hand every slot to the two or three strongest
    chassis, and the strongest switch is not the only thing a slot is bought for: a design whose
    parent the novelty gate rejects scores nothing at all, so the chassis predicted most novel
    has to be on the plate even though its switch is the weakest measured here. Taking one round
    at a time gives every admissible parent a slot before any parent gets a second, and the
    rounds after that still go in switch order.
    """
    chosen: list[Candidate] = []
    remaining = group_by(ranked, lambda candidate: candidate.parent)
    while len(chosen) < slots:
        before = len(chosen)
        for parent in sorted(remaining, key=lambda p: -remaining[p][0].switch):
            queue = remaining[parent]
            at = next((i for i, c in enumerate(queue)
                       if not any(c.too_similar_to(taken) for taken in chosen)), None)
            if at is not None:
                chosen.append(queue.pop(at))
        remaining = {parent: rest for parent, rest in remaining.items() if rest}
        if len(chosen) == before:
            break
    return chosen[:slots]


def report_parents(parents: Sequence[str], novelty: Mapping[str, Novelty],
                   passes: Mapping[str, bool], unpaired: Mapping[str, float]) -> None:
    typer.echo(f"{'parent':34s} {'novelty':>7s} {'triage':>6s} {'unpairedC':>9s}  note")
    for parent in parents:
        typer.echo(f"{parent:34s} {novelty[parent].score:7d} "
                   f"{'pass' if passes.get(parent) else 'FAIL':>6s} "
                   f"{unpaired.get(parent, 0.0):9.0f}  {novelty[parent].note}")


def report_chosen(rows: Sequence[Mapping[str, object]]) -> None:
    typer.echo(f"\n{'rank':>4s} {'name':44s} {'logU':>5s} {'fold':>6s} {'sites':>5s} "
               f"{'epi':>12s} {'nov':>3s}")
    for row in rows:
        typer.echo(f"{row['rank']:4d} {str(row['name']):44s} "
                   f"{row['switch_log_units']:5.2f} {row['fold_weaker_at_ph6']:6.1f} "
                   f"{row['coupled_sites']:5d} {str(row['epitope']):>12s} "
                   f"{row['novelty_score']:3d}")


@app.command()
def main(
    switch: Annotated[Path, typer.Option(help="per design/position agreement CSV")],
    variants: Annotated[Path, typer.Option(help="name,sequence CSV the variants were built as")],
    pooled: Annotated[Path, typer.Option(help="pooled per-design metrics")],
    triage: Annotated[Path, typer.Option(help="per-design triage verdicts")],
    sequence_hits: Annotated[Path, typer.Option(help="MMseqs2 m8 over the reference corpus")],
    fold_hits: Annotated[Path, typer.Option(help="Foldseek m8 over CATH50")],
    conservation: Annotated[Path, typer.Option(help="measured epitope hotspots JSON")],
    out: Annotated[Path, typer.Option(help="the chosen designs")],
    rejected: Annotated[Path, typer.Option(help="every candidate with its reason")],
    # Proteinbase's cap is 20 per submission, stated once with its provenance in
    # `tnf_binder.submission.MAX_ROWS`; this cannot import it, so it must not drift from it.
    slots: Annotated[int, typer.Option(help="submission cap for the track")] = 20,
    novelty_floor: Annotated[int, typer.Option(help="lowest Proteinbase band to submit")] = 3,
    allow_unpaired_cys: Annotated[bool, typer.Option(
        help="admit chassis carrying an unpaired cysteine")] = False,
    thresholds: Annotated[Path | None, typer.Option(
        help="the verifier's <out_summary>.thresholds.json, which states the floor it used")] = None,
    floor: Annotated[float | None, typer.Option(
        "--coupling-floor",
        help="override the linkage floor instead of reading it from --thresholds")] = None,
) -> None:
    coupling = coupling_floor(thresholds, floor)
    sequence_m8 = read_m8(sequence_hits,
                          ("query", "target", "pident", "alnlen", "qcov", "evalue", "bits"))
    fold_m8 = read_m8(fold_hits, ("query", "target", "alntmscore", "qtmscore", "ttmscore",
                                  "qcov", "tcov", "evalue"))
    # Both the epitope of each design and its mouse identity are measurements already on disk:
    # `collect_designs` writes the epitope per design, and `prepare_target` writes each patch's
    # measured conservation. Re-deriving either by prefix-matching design names would be a
    # second, divergeable copy of a number the funnel is reported against.
    metrics = {row["name"]: row for row in read_rows(pooled)}
    passes = {row["name"]: row["passes"] == "yes" for row in read_rows(triage)}
    unpaired = {name: number(row.get("unpaired_cys")) for name, row in metrics.items()}
    conserved = {name: float(patch["conservation"])
                 for name, patch in json.loads(conservation.read_text())["epitopes"].items()}

    pool = candidates(read_rows(switch), read_rows(variants), coupling)
    parents = sorted({candidate.parent for candidate in pool})
    novelty = {parent: novelty_of(parent, sequence_m8, fold_m8) for parent in parents}

    typer.echo(f"{len(pool)} variants over {len(parents)} parent designs, "
               f"coupling floor {coupling}\n")
    report_parents(parents, novelty, passes, unpaired)

    reasons = {c.name: admissible(c, novelty, passes, novelty_floor, allow_unpaired_cys, unpaired)
               for c in pool}
    eligible = [c for c in pool if not reasons[c.name]]
    ranked = sorted(eligible, key=lambda c: (
        -c.switch,
        -conserved.get(metrics.get(c.parent, {}).get("epitope", ""), 0.0),
        -number(metrics.get(c.parent, {}).get("iptm")),
    ))
    chosen = choose(ranked, slots)

    rows = [c.row(metrics.get(c.parent, {}), novelty[c.parent], conserved) | {"rank": position}
            for position, c in enumerate(chosen, 1)]
    write_rows(out, rows, OUT_COLUMNS)

    taken = {c.name for c in chosen}
    write_rows(rejected, [
        {"name": c.name, "parent": c.parent, "switch_log_units": round(c.switch, 2),
         "coupled_sites": c.coupled, "chosen": c.name in taken,
         "reason": reasons[c.name] or ("" if c.name in taken else "not in the top slots")}
        for c in sorted(pool, key=lambda c: (-c.switch, c.name))
    ], REJECTED_COLUMNS)

    typer.echo(f"\n{len(eligible)} of {len(pool)} variants are admissible; {len(chosen)} chosen")
    report_chosen(rows)
    typer.echo(f"\nwrote {out} and {rejected}")


if __name__ == "__main__":
    app()
