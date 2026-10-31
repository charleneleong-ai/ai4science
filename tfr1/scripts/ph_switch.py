"""Is each histidine's titration actually coupled to binding? Two methods, one pass.

A histidine whose pKa sits in the 6.0-7.6 window is *necessary* for a pH switch and nowhere
near sufficient. If it titrates identically whether the binder is bound or free, the
protonation equilibrium is the same on both sides and cancels: the site contributes nothing
to pH-dependent binding. What governs the pH dependence of Kd is the linkage between the two
states,

    delta = pKa(complex) - pKa(free)

which must be meaningfully **negative** for the switch this challenge asks for -- protonation
favoured in the free state, so dropping the pH drives dissociation. A positive delta binds
*tighter* at low pH, which is backwards.

The TNF-alpha campaign learned this the expensive way: it ranked an entry on six histidines
in the window, then measured the linkage and found its two best-determined positions
uncoupled on all 25 samples, all six sitting 9-21 A from the target with none in contact.
Four things it did that this does differently:

- **It measured the linkage with PROPKA alone.** PROPKA is an empirical scheme; if its
  systematic error points the same way at every histidine on one interface, the verdict
  inherits that error and no amount of extra structures reveals it. Here both PROPKA and
  pKAI -- a graph neural network trained on experimental values, sharing no machinery with
  it -- compute delta, and disagreement is reported rather than averaged away. Validating
  this script against TNF's own structures found PROPKA calling a non-contacting histidine
  coupled (delta -0.41) where pKAI put it at -0.11, which is the case a single method cannot
  see.
- **It ranked on a count of in-window histidines.** What a switch is worth is a free energy,
  so `release_log_units` converts each site's linkage into the log-fold weakening of Kd
  between the two assay pH values, and designs are ranked on the sum. One site buys little:
  the measured positive case is 0.2-0.4 log units and the physical ceiling is about 1.4, so
  "no detectable binding at pH 6.0" needs three or four coupled sites.
- **It judged each position on a mean.** Across independent structures of one design the
  geometry varies enough to carry the mean on its own: TNF's His64 sat at 13.3, 6.9 and 5.3 A
  in three models with deltas of -0.00, -0.07 and -1.16, so the mean says "coupled" on the
  strength of one outlier. Positions here are judged on the median, and a spread wider than
  the coupling floor is reported as `inconsistent` rather than resolved either way.
- **The target's sidechains were never checked.** BoltzGen reproduces its input target
  faithfully, truncation included: the GLP1R target it was given was missing the carboxylate
  of 9 of 48 residues (ASP at 5 heavy atoms is N/CA/C/O/CB -- no OD1/OD2). A missing
  carboxylate is invisible to both pKa methods and drives delta toward zero, manufacturing
  exactly the uncoupled verdict TNF reported. `prepare_target.py` now refuses to *write*
  such a target; this refuses to measure one, because a later generator can drop atoms its
  input had.

The apo state is the complex **minus the target chain**, holding the bound conformation.
Predicting the free binder separately would conflate conformational relaxation with the
electrostatic linkage, which is not what delta means.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import math
import shutil
import statistics
import tempfile
from collections.abc import Callable, Collection, Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Protocol, TypeVar

import gemmi
import typer
from csvio import read_designs, write_rows
from target_sites import (
    CARBOXYLATE_ATOMS,
    one_letter,
    placeholder_coordinates,
    truncated_carboxylates,
)

# The window the previous campaign gated on. Kept because recording a rejected metric is what
# makes its rejection reproducible, not because anything here ranks on it -- see
# `release_log_units` for why gating on the bound-state window discards working switches.
WINDOW = (6.0, 7.6)

VERDICT_RELEASE = "release at low pH"
VERDICT_WRONG = "wrong direction"
VERDICT_UNCOUPLED = "uncoupled"
VERDICT_INCONSISTENT = "inconsistent across structures"
VERDICT_UNREPLICATED = "single structure - not replicated"

RANK_KEY = "total_release_consistent"

K = TypeVar("K")
T = TypeVar("T")


def group_by(rows: Iterable[T], key: Callable[[T], K]) -> dict[K, list[T]]:
    grouped: dict[K, list[T]] = {}
    for row in rows:
        grouped.setdefault(key(row), []).append(row)
    return grouped


def release_log_units(pka_free: float, pka_bound: float, bind: float, release: float) -> float:
    """How much weaker Kd gets between the binding pH and the release pH, in log10.

    For a single titratable site the linkage relation gives

        log Kd(pH) = log Kd(high pH) + log10[(1 + 10^(pKa_free - pH))
                                             / (1 + 10^(pKa_bound - pH))]

    so what matters is the difference of that term between the two assay pH values. Positive
    means weaker binding at the release pH, which is the switch being asked for.

    This is the statistic to rank on, and it is *not* the one the TNF-alpha campaign gated on.
    That gate asked whether the **complex** pKa fell inside 6.0-7.6, and validating this
    script on TNF's own data rejects that gate: its single cleanly-coupled histidine had a
    complex pKa of 3.4-5.7, well below the window, while its apo and complex midpoints
    straddled the release pH exactly as a working switch must. Gating on the bound-state
    window would have thrown away the one site that worked.
    """
    def term(pka: float, ph: float) -> float:
        return math.log10(1.0 + 10.0 ** (pka - ph))

    return ((term(pka_free, release) - term(pka_bound, release))
            - (term(pka_free, bind) - term(pka_bound, bind)))


@dataclass(frozen=True)
class Limits:
    """Every number that decides a verdict, resolved from the CLI and written beside the
    results. `target_sites.py` does the same, for the same reason: an earlier revision of this
    campaign hardcoded its thresholds and shipped two epitopes its own gate had rejected."""

    coupling_floor: float
    ph_bind: float
    ph_release: float

    def as_dict(self) -> dict[str, object]:
        return {
            "coupling_floor": self.coupling_floor,
            "ph_bind": self.ph_bind,
            "ph_release": self.ph_release,
            "in_window_reported_not_ranked": list(WINDOW),
            "rank_key": RANK_KEY,
        }


@dataclass(frozen=True)
class Reading:
    """One histidine, one structure, one method."""

    design: str
    source: str
    method: str
    position: int
    pka_complex: float
    pka_apo: float
    target_distance: float
    nearest_carboxylate: str
    limits: Limits

    @property
    def delta(self) -> float:
        return self.pka_complex - self.pka_apo

    @property
    def coupled(self) -> bool:
        return self.delta <= -self.limits.coupling_floor

    @property
    def in_window(self) -> bool:
        return WINDOW[0] <= self.pka_complex <= WINDOW[1]

    @property
    def release(self) -> float:
        return release_log_units(self.pka_apo, self.pka_complex,
                                 self.limits.ph_bind, self.limits.ph_release)

    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "source": self.source,
            "method": self.method,
            "position": self.position,
            "pka_complex": round(self.pka_complex, 2),
            "pka_apo": round(self.pka_apo, 2),
            "delta": round(self.delta, 2),
            "release_log_units": round(self.release, 3),
            "coupled": self.coupled,
            "in_window": self.in_window,
            "target_distance": round(self.target_distance, 1),
            "nearest_carboxylate": self.nearest_carboxylate,
        }


@dataclass(frozen=True)
class Summary:
    """One histidine under one method, across every structure of its design."""

    design: str
    position: int
    method: str
    n_structures: int
    delta_median: float
    delta_spread: float
    release_median: float
    coupled_count: int
    in_window_count: int
    target_distance: float
    nearest_carboxylate: str
    verdict: str

    @property
    def releasing(self) -> bool:
        return self.verdict == VERDICT_RELEASE

    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "position": self.position,
            "method": self.method,
            "n_structures": self.n_structures,
            "delta_median": round(self.delta_median, 2),
            "delta_spread": round(self.delta_spread, 2),
            "release_median": round(self.release_median, 3),
            "coupled": f"{self.coupled_count}/{self.n_structures}",
            "in_window": f"{self.in_window_count}/{self.n_structures}",
            "target_distance": round(self.target_distance, 1),
            "nearest_carboxylate": self.nearest_carboxylate,
            "verdict": self.verdict,
        }


class Method(Protocol):
    name: str

    def pkas(self, pdb: Path) -> dict[tuple[str, int], float]: ...


class Propka:
    """Empirical: desolvation, hydrogen bonding and Coulomb terms with fitted coefficients."""

    name = "propka"

    def pkas(self, pdb: Path) -> dict[tuple[str, int], float]:
        import propka.run  # local: heavy, and the pure functions above must import without it

        # PROPKA narrates every imperfect sidechain it meets, on both streams.
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            molecule = propka.run.single(str(pdb), optargs=["--quiet"], write_pka=False)
        conformation = molecule.conformations[molecule.conformation_names[0]]
        return {
            (group.atom.chain_id, int(group.atom.res_num)): float(group.pka_value)
            for group in conformation.groups
            if group.residue_type == "HIS"
        }


class PKai:
    """A graph neural network trained on experimental and PHMD values. Shares no machinery
    with PROPKA, which is the whole point of asking it."""

    name = "pkai"

    def pkas(self, pdb: Path) -> dict[tuple[str, int], float]:
        from pkai.pKAI import pKAI  # local: pulls in torch

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            raw = pKAI(str(pdb), device="cpu")
        return {
            (chain, int(number)): float(value)
            for chain, number, name, value in raw
            if name == "HIS"
        }


class TargetUnmeasurable(Exception):
    """A structure that cannot honestly be measured -- distinct from a CLI usage error."""


def chain_sequence(chain: gemmi.Chain) -> str:
    """The chain's one-letter sequence under the convention the rest of this package uses.

    `target_sites.one_letter` uppercases and drops anything that is not an amino acid, where
    `gemmi.one_letter_code` lowercases non-standard residues. The recorded sequences this is
    compared against are uppercased, so the conventions have to match or a perfectly good
    structure fails as "wrong chain".
    """
    return "".join(one_letter(residue) for residue in chain)


def target_atoms(chains: Collection[gemmi.Chain]) -> list[tuple[gemmi.Position, str, bool]]:
    """Every target atom once per structure, tagged with its residue and whether that residue
    carries a carboxylate. Built once and shared by every histidine, because walking the gemmi
    hierarchy per histidine is the expensive part of the geometry -- flattening it first
    measured 11-18x faster on this campaign's own targets."""
    return [
        (atom.pos, f"{residue.name}{residue.seqid.num}", residue.name in CARBOXYLATE_ATOMS)
        for chain in chains for residue in chain for atom in residue
    ]


def geometry_of(residue: gemmi.Residue,
                atoms: Sequence[tuple[gemmi.Position, str, bool]]) -> tuple[float, str]:
    """Distance to the nearest target atom and to the nearest receptor carboxylate, in one
    pass. Both are the mechanistic covariates that explain a coupling result, so both come
    from the same sweep."""
    probes = [atom.pos for atom in residue]
    nearest = float("inf")
    acid_name, acid_distance = "none", float("inf")
    for position, label, acidic in atoms:
        distance = min(probe.dist(position) for probe in probes)
        if distance < nearest:
            nearest = distance
        if acidic and distance < acid_distance:
            acid_name, acid_distance = label, distance
    label = "none" if acid_distance == float("inf") else f"{acid_name}@{acid_distance:.1f}A"
    return nearest, label


@dataclass(frozen=True)
class Complex:
    """A loaded complex that is admissible to measure.

    Everything that can make a structure unmeasurable is checked in `load_complex`, so the
    existence of this object means the binder is identified and the target is complete. The
    checks deliberately do not live in `measure`, which is a generator: a gate that fires only
    when the caller happens to iterate is not a gate.
    """

    design: str
    source: str
    binder: gemmi.Chain
    targets: tuple[gemmi.Chain, ...]
    structure: gemmi.Structure

    @property
    def histidines(self) -> list[gemmi.Residue]:
        return [residue for residue in self.binder if residue.name == "HIS"]

    def write_pair(self, work: Path) -> tuple[Path, Path]:
        complex_path = work / f"{self.design}_{self.source}_complex.pdb"
        self.structure.write_pdb(str(complex_path))

        apo = gemmi.Structure()
        apo.add_model(gemmi.Model("1"))
        apo[0].add_chain(self.binder)
        apo.setup_entities()
        apo_path = work / f"{self.design}_{self.source}_apo.pdb"
        apo.write_pdb(str(apo_path))
        return complex_path, apo_path


def plain_copy(path: Path, work: Path) -> Path:
    """Structures are gzipped to reclaim disk; both pKa tools need them plain."""
    if path.suffix != ".gz":
        return path
    out = work / path.stem
    with gzip.open(path, "rb") as source, out.open("wb") as destination:
        shutil.copyfileobj(source, destination)
    return out


def load_complex(path: Path, design: str, sequence: str, work: Path) -> Complex:
    structure = gemmi.read_structure(str(plain_copy(path, work)))
    structure.setup_entities()
    structure.remove_ligands_and_waters()
    model = structure[0]

    binders = [chain for chain in model
               if len(chain) == len(sequence) and chain_sequence(chain) == sequence]
    if not binders:
        shapes = {chain.name: f"{len(chain)} res" for chain in model}
        raise TargetUnmeasurable(
            f"{path.name}: no chain matches the recorded sequence for {design} "
            f"({len(sequence)} res); chains are {shapes}. Refusing to guess which chain is "
            f"the binder -- measuring the wrong one would poison every delta."
        )
    if len(binders) > 1:
        raise TargetUnmeasurable(
            f"{path.name}: {len(binders)} chains carry {design}'s sequence, so the binder is "
            f"ambiguous and the duplicate would be scored as receptor."
        )

    binder = binders[0]
    targets = tuple(chain for chain in model if chain.name != binder.name)
    if not targets:
        raise TargetUnmeasurable(f"{path.name}: no target chain, so there is no linkage")

    residues = [residue for chain in targets for residue in chain]
    if missing := truncated_carboxylates(residues):
        raise TargetUnmeasurable(
            f"{path.name}: {len(missing)} target carboxylates are missing their charged atoms "
            f"({', '.join(missing[:5])}). Both pKa methods would score them as absent and "
            f"drive delta toward zero, manufacturing an uncoupled verdict."
        )

    everything = [residue for chain in model for residue in chain]
    if collapsed := placeholder_coordinates(everything):
        raise TargetUnmeasurable(
            f"{path.name}: {len(collapsed)} residues have their sidechain atoms parked at the "
            f"origin rather than placed ({', '.join(collapsed[:5])}). The atoms are all present "
            f"by name, so nothing downstream would notice; both pKa methods would score a "
            f"sidechain collapsed onto one point and return a reading rather than fail."
        )

    source = path.name.removesuffix(".gz").removesuffix(".cif")
    return Complex(design, source, binder, targets, structure)


def index_structures(root: Path, designs: Iterable[str]) -> dict[str, list[Path]]:
    """One walk of the tree, not one per design.

    Generators disagree about where the design name sits in the filename: BoltzGen writes
    `rank04_conserved_17.cif`, Boltz-2 writes `<design>_model_0.cif.gz`. Both are accepted by
    testing for the name delimited by underscores, so `conserved_1` cannot claim
    `rank01_conserved_11`, and the longest match wins so the result is order-independent.
    """
    found: dict[str, list[Path]] = {name: [] for name in designs}
    for path in sorted(root.rglob("*.cif*")):
        if path.suffixes[-2:] not in ([".cif"], [".cif", ".gz"]):
            continue
        stem = path.name.removesuffix(".gz").removesuffix(".cif")
        padded = f"_{stem}_"
        if matches := [name for name in found if f"_{name}_" in padded]:
            found[max(matches, key=len)].append(path)
    return found


def measure(entry: Complex, methods: Iterable[Method], limits: Limits,
            work: Path) -> Iterator[Reading]:
    histidines = entry.histidines
    if not histidines:
        return

    complex_path, apo_path = entry.write_pair(work)
    atoms = target_atoms(entry.targets)
    geometry = {residue.seqid.num: geometry_of(residue, atoms) for residue in histidines}
    for method in methods:
        bound, free = method.pkas(complex_path), method.pkas(apo_path)
        for residue in histidines:
            position = residue.seqid.num
            pair = (entry.binder.name, position)
            if pair not in bound or pair not in free:
                continue
            distance, carboxylate = geometry[position]
            yield Reading(
                design=entry.design,
                source=entry.source,
                method=method.name,
                position=position,
                pka_complex=bound[pair],
                pka_apo=free[pair],
                target_distance=distance,
                nearest_carboxylate=carboxylate,
                limits=limits,
            )


def verdict(group: Sequence[Reading], limits: Limits) -> str:
    """Judge a position on whether its structures *agree*, never on their mean.

    One outlier structure is enough to drag a mean across the floor -- TNF's His64 moved
    between 13.3 and 5.3 A across three models, with deltas of -0.00, -0.07 and -1.16, so its
    mean says "coupled" on the strength of one model out of three. So a verdict is only
    returned when every structure gives the same answer, and a position where some structures
    clear the floor and others do not is `inconsistent`.

    The test is deliberately on the verdict and not on the magnitude spread: -1.48, -1.65 and
    -2.64 vary by more than the floor while agreeing unanimously that the site is coupled,
    which is a real switch with geometry-dependent strength, not a disagreement. A position
    seen in a single structure is not replicated and says nothing about consistency either
    way. `delta_spread` is still reported, because the magnitude variation is informative --
    it is just not what decides the verdict.
    """
    if len(group) < 2:
        return VERDICT_UNREPLICATED
    coupled = sum(1 for reading in group if reading.coupled)
    wrong = sum(1 for reading in group if reading.delta >= limits.coupling_floor)
    if coupled == len(group):
        return VERDICT_RELEASE
    if wrong == len(group):
        return VERDICT_WRONG
    return VERDICT_INCONSISTENT if coupled or wrong else VERDICT_UNCOUPLED


def summarise(readings: Iterable[Reading], limits: Limits) -> list[Summary]:
    grouped = group_by(readings, lambda r: (r.design, r.position, r.method))
    out: list[Summary] = []
    for (design, position, method), group in sorted(grouped.items()):
        deltas = [reading.delta for reading in group]
        out.append(Summary(
            design=design,
            position=position,
            method=method,
            n_structures=len(group),
            delta_median=statistics.median(deltas),
            delta_spread=max(deltas) - min(deltas),
            release_median=statistics.median(reading.release for reading in group),
            coupled_count=sum(1 for reading in group if reading.coupled),
            in_window_count=sum(1 for reading in group if reading.in_window),
            target_distance=statistics.mean(reading.target_distance for reading in group),
            nearest_carboxylate=group[0].nearest_carboxylate,
            verdict=verdict(group, limits),
        ))
    return out


def agreement(summary: Iterable[Summary]) -> list[dict[str, object]]:
    """Where the two methods disagree about a position, say so rather than average."""
    grouped = group_by(summary, lambda s: (s.design, s.position))
    out: list[dict[str, object]] = []
    for (design, position), rows in sorted(grouped.items()):
        ordered = sorted(rows, key=lambda row: row.method)
        first = ordered[0]
        out.append({
            "design": design,
            "position": position,
            **{f"delta_{row.method}": round(row.delta_median, 2) for row in ordered},
            **{f"release_{row.method}": round(row.release_median, 3) for row in ordered},
            "verdicts": " | ".join(f"{row.method}:{row.verdict}" for row in ordered),
            "methods_agree": len({row.verdict for row in ordered}) == 1,
            "both_releasing": len(ordered) > 1 and all(row.releasing for row in ordered),
            "target_distance": round(first.target_distance, 1),
            "nearest_carboxylate": first.nearest_carboxylate,
        })
    return out


def per_design(summary: Iterable[Summary]) -> list[dict[str, object]]:
    """Add each design's sites up, because in the independent-site approximation the log-unit
    contributions add. Two totals are reported: every site, and only the sites whose verdict
    survived the consistency check. The second is what designs are ranked on, so a total
    carried by one outlier structure cannot win."""
    out: list[dict[str, object]] = []
    for (design, method), rows in group_by(summary, lambda s: (s.design, s.method)).items():
        consistent = sum(row.release_median for row in rows if row.releasing)
        out.append({
            "design": design,
            "method": method,
            "n_his": len(rows),
            "n_releasing": sum(1 for row in rows if row.releasing),
            "n_inconsistent": sum(1 for row in rows if row.verdict == VERDICT_INCONSISTENT),
            "n_contacting": sum(1 for row in rows if row.target_distance <= 5.0),
            "total_release_all_sites": round(sum(row.release_median for row in rows), 2),
            RANK_KEY: round(consistent, 2),
            "fold_weaker_at_ph6": round(10.0 ** consistent, 1),
        })
    return sorted(out, key=lambda row: float(row[RANK_KEY]), reverse=True)


def report_positions(compared: Sequence[dict[str, object]]) -> None:
    typer.echo(f"\n{'design':28s} {'His':>4s} {'dprop':>6s} {'dpkai':>6s} {'dist':>5s} "
               f"{'agree':>5s}  verdicts")
    for row in compared:
        propka = row.get("delta_propka", float("nan"))
        pkai = row.get("delta_pkai", float("nan"))
        typer.echo(f"{row['design']:28s} {row['position']:4d} {propka:6.2f} {pkai:6.2f} "
                   f"{row['target_distance']:5.1f} "
                   f"{'yes' if row['methods_agree'] else 'NO':>5s}  {row['verdicts']}")


def report_designs(ranked: Sequence[dict[str, object]]) -> None:
    typer.echo(f"\n{'design':28s} {'method':7s} {'His':>4s} {'rel':>4s} {'incon':>6s} "
               f"{'contact':>8s} {'logKd':>6s} {'fold':>6s}")
    for row in ranked:
        typer.echo(f"{row['design']:28s} {row['method']:7s} {row['n_his']:4d} "
                   f"{row['n_releasing']:4d} {row['n_inconsistent']:6d} "
                   f"{row['n_contacting']:8d} {row[RANK_KEY]:6.2f} "
                   f"{row['fold_weaker_at_ph6']:6.1f}")
    typer.echo("\nA switch that reads as 'no detectable binding at pH 6.0' needs roughly "
               "2 log units, so three or four coupled sites. One is not enough.")


def collect(sequences: dict[str, str], structures: Path, methods: Sequence[Method],
            limits: Limits) -> tuple[list[Reading], dict[str, int]]:
    index = index_structures(structures, sequences)
    readings: list[Reading] = []
    tally = {"no_structure": 0, "no_histidine": 0, "measured": 0, "unmeasurable": 0}
    for design, sequence in sequences.items():
        found = index[design]
        if not found:
            tally["no_structure"] += 1
            continue
        new: list[Reading] = []
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            for path in found:
                try:
                    entry = load_complex(path, design, sequence, work)
                except TargetUnmeasurable as exc:
                    typer.echo(f"  SKIPPED {exc}")
                    tally["unmeasurable"] += 1
                    continue
                new.extend(measure(entry, methods, limits, work))
        if new:
            readings.extend(new)
            tally["measured"] += 1
            typer.echo(f"  {design:28s} {len(found)} structures, {len(new)} readings")
        else:
            tally["no_histidine"] += 1
    return readings, tally


app = typer.Typer(add_completion=False)


@app.command()
def main(
    designs: Annotated[Path, typer.Option(help="pooled CSV from collect_designs.py")],
    structures: Annotated[Path, typer.Option(help="root to search for complex CIFs")],
    out_rows: Annotated[Path, typer.Option(help="one row per design/structure/method/His")],
    out_summary: Annotated[Path, typer.Option(help="per design/His/method summary")],
    out_agreement: Annotated[Path, typer.Option(help="per design/His cross-method comparison")],
    out_designs: Annotated[Path, typer.Option(help="per design/method summed switch strength")],
    only: Annotated[str, typer.Option(help="comma-separated design names to restrict to")] = "",
    coupling_floor: Annotated[float, typer.Option(
        help="how far a pKa must drop on binding before coupling is called real")] = 0.3,
    ph_bind: Annotated[float, typer.Option(help="the pH binding is measured at")] = 7.4,
    ph_release: Annotated[float, typer.Option(help="the pH release is measured at")] = 6.0,
) -> None:
    limits = Limits(coupling_floor=coupling_floor, ph_bind=ph_bind, ph_release=ph_release)
    wanted = {name.strip() for name in only.split(",") if name.strip()} or None
    sequences = read_designs(designs, wanted)
    methods: list[Method] = [Propka(), PKai()]

    readings, tally = collect(sequences, structures, methods, limits)
    typer.echo(f"\n{len(sequences)} designs: {tally['no_structure']} with no structure on disk, "
               f"{tally['no_histidine']} carrying no histidine, {tally['measured']} measured "
               f"({tally['unmeasurable']} structures skipped as unmeasurable)")
    if not readings:
        typer.echo("no histidine in any binder chain -- nothing to couple. The pH switch has "
                   "to be engineered in before this stage can say anything.")
        return

    summary = summarise(readings, limits)
    compared = agreement(summary)
    ranked = per_design(summary)

    write_rows(out_rows, [reading.row() for reading in readings])
    write_rows(out_summary, [row.row() for row in summary])
    write_rows(out_agreement, compared)
    write_rows(out_designs, ranked)
    sidecar = out_summary.with_suffix(".thresholds.json")
    sidecar.write_text(json.dumps(limits.as_dict(), indent=2) + "\n")

    report_positions(compared)
    coupled = [row for row in compared if row["both_releasing"]]
    typer.echo(f"\n{len(coupled)} of {len(compared)} histidine positions are coupled in the "
               f"releasing direction with both methods agreeing")
    report_designs(ranked)
    typer.echo(f"wrote {out_rows}, {out_summary}, {out_agreement}, {out_designs} and {sidecar}")


if __name__ == "__main__":
    app()
