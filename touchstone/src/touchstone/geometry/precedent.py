"""Coordination-motif precedent check — has this metal environment been seen before?

A cheap, strong prior: a metal–donor motif with many precedents is well-precedented (trust);
an unprecedented one (no hits) sits off the known manifold and is deferred for deeper scrutiny.
As a *cheap* gate (a table lookup, no GPU) it belongs at the front of `cascade`, so the
expensive tiers only run on designs with structural precedent.

The default searcher is **open**: `metalpdb_precedent_search` reads a bundled motif→count table mined
from MetalPDB (scripts/build_metalpdb_precedents.py) — the licence-free, metalloprotein-specific
analog of a CSD search. A CSD-CrossMiner searcher (`crossminer_search`, licence-gated) is
provided for CSD users. The searcher is pluggable (inject one for tests or a different source).
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable

from ..core import BinderDesign, CoordinationSite, Verdict, element_symbol

PRECEDENTS_DATA = Path(__file__).parent.parent / "data" / "metalpdb_precedents.json"


@dataclass
class PrecedentHits:
    """Precedent support for one coordination motif."""

    nhits: int  # structures matching the metal–donor environment
    motif: str  # e.g. "Ni-N3O2"


def motif_label(site: CoordinationSite) -> str:
    """A compact label for a coordination environment, e.g. 'Ni-N3O2'."""
    return f"{element_symbol(site.metal)}-{donor_set(site)}"


@lru_cache(maxsize=1)
def precedents_table() -> dict[str, int]:
    """The bundled MetalPDB motif→count table (built by scripts/build_metalpdb_precedents.py)."""
    return json.loads(PRECEDENTS_DATA.read_text())


def metalpdb_precedent_search(site: CoordinationSite) -> PrecedentHits:
    """Open precedent via the bundled MetalPDB motif→count table — licence-free, the
    metalloprotein analog of a CSD-CrossMiner search. Raises if the table isn't bundled."""
    label = motif_label(site)
    return PrecedentHits(nhits=int(precedents_table().get(label, 0)), motif=label)


def crossminer_search(site: CoordinationSite) -> PrecedentHits:
    """Search CSD-CrossMiner for the coordination motif. Requires `ccdc` + a licence;
    the pharmacophore/connectivity query should be confirmed against the installed API."""
    from ccdc.search import SubstructureSearch  # noqa: F401 — licence probe

    raise NotImplementedError(
        f"wire the CrossMiner query for {motif_label(site)} against the licensed CSD API"
    )


MIN_METAL_SITES = 30  # a metal with fewer precedented sites can't support an enrichment (Au3+ has 4)
MIN_MOTIF_HITS = 3  # ... and needs this many sites *with this donor set* before we call it the owner
SMOOTHING = 5.0  # pseudo-observations pulling a thin metal toward the background rate, not past it


@lru_cache(maxsize=1)
def metal_totals() -> dict[str, int]:
    """Total precedented sites per metal element — the denominator for enrichment.

    Raw hit counts are useless as a metal comparison: Zn has ~10× more sites in the PDB than Ni, so
    it would look 'preferred' for every motif on abundance alone. Normalising by each metal's total
    presence is what turns the counts into chemistry."""
    totals: Counter[str] = Counter()
    for motif, n in precedents_table().items():
        totals[motif.split("-", 1)[0]] += n
    return dict(totals)


@dataclass
class MotifEnrichment:
    """How characteristic a donor set is of each metal, relative to that metal's whole PDB presence."""

    donors: str  # the donor set alone, e.g. "N2S2"
    target: str
    enrichment: dict[str, float]  # metal → smoothed fraction of that metal's sites using this set
    hits: dict[str, int]  # raw counts — a metal with < MIN_MOTIF_HITS cannot own the donor set
    background: float  # the pooled rate across the panel: what "unremarkable" looks like

    @property
    def eligible(self) -> list[str]:
        """Metals with enough observed sites on this donor set to be called its owner. Without this,
        a metal that has *never* been seen on the donor set can still top the ranking on smoothing
        alone — the enrichment is a rate, and a rate estimated from zero observations is not
        evidence of anything."""
        return [m for m, n in self.hits.items() if n >= MIN_MOTIF_HITS]

    @property
    def preferred(self) -> str | None:
        return max(self.eligible, key=lambda m: self.enrichment[m], default=None)

    @property
    def ratio(self) -> float:
        """Target enrichment ÷ the best competitor's. >1 ⇒ the donor set is more characteristic of
        the target metal than of anything else in the panel. Competitors are smoothed, never zero,
        so this is always finite — but it needs at least one competitor to mean anything."""
        competitors = [e for m, e in self.enrichment.items() if m != self.target]
        best = max(competitors, default=0.0)
        return self.enrichment[self.target] / best if best else 1.0


def donor_set(site: CoordinationSite) -> str:
    """The donor composition alone, e.g. 'N2S2' — the metal-free half of a motif label."""
    return "".join(f"{el}{n}" for el, n in sorted(Counter(site.ligand_elems).items()))


def motif_enrichment(site: CoordinationSite, metals: tuple[str, ...]) -> MotifEnrichment:
    """Enrichment of this site's donor set across a panel of metals, from the bundled MetalPDB table.

    Smoothed toward the **background rate** for this donor set (its pooled frequency over the panel),
    not toward a constant. Add-one smoothing would be wrong here: each metal's denominator is its own
    total, and those span 60× (Pt 53 → Zn 3105), so a constant pseudo-count is a *metal-dependent
    floor* — it hands a rare metal with zero hits a higher rate than an abundant metal with real
    ones. Shrinking toward the pooled rate instead means a thin metal regresses to unremarkable,
    which is the honest answer when you have no observations."""
    table, totals = precedents_table(), metal_totals()
    panel = tuple(dict.fromkeys((site.metal, *metals)))
    donors = donor_set(site)

    hits = {m: int(table.get(f"{element_symbol(m)}-{donors}", 0)) for m in panel}
    sites = {m: totals.get(element_symbol(m), 0) for m in panel}
    pooled = sum(sites.values())
    background = sum(hits.values()) / pooled if pooled else 0.0

    enrichment = {
        m: (hits[m] + SMOOTHING * background) / (sites[m] + SMOOTHING) for m in panel
    }
    return MotifEnrichment(donors, site.metal, enrichment, hits, background)


class MotifSelectivityVerifier:
    """Does this donor set actually belong to the target metal, in real proteins?

    The metal discrimination the geometry tiers can't make and the MLIP tier can't be trusted to make
    (see docs/experiments/2026-07-13-mlip-cannot-rank-metals.md) — recovered instead from *observed
    occupancy*: of all the sites a metal occupies in the PDB, what fraction use this donor set?
    Normalising is the whole trick (see `metal_totals`).

    Where the evidence is there, it recovers HSAB without computing any physics: O6 is most
    characteristic of Mn²⁺ (hard) and least of Cu²⁺; N2S2 (His₂/Cys/Met) of Cu²⁺ (soft, type-1
    blue copper); S4 of Fe (iron–sulfur clusters). Where it isn't, it defers.

    **Scope — this is an occupancy prior, not a binding free energy.** It says a donor set is
    *characteristic* of a metal in biology, which is confounded by cellular abundance and by what
    evolution needed. It cannot tell you that Cu²⁺ out-competes Ni²⁺ in a mixed leachate — that is a
    thermodynamic question, and answering it needs ligand-field physics no MLIP we have carries.
    Useful for *designing a site that looks like a real Cu site*; not a substitute for ΔΔG."""

    def __init__(
        self,
        metals: tuple[str, ...] = ("Ni2+", "Cu2+", "Co2+"),
        *,
        trust_ratio: float = 1.25,  # target must be this many × more enriched than the best competitor
        min_metal_sites: int = MIN_METAL_SITES,
    ) -> None:
        self.metals = metals
        self.trust_ratio = trust_ratio
        self.min_metal_sites = min_metal_sites

    def panel(self, design: BinderDesign) -> list[str]:
        """The target plus its competitors, one per *element*: the table keys on the element, so
        Fe2+ and Fe3+ resolve to the same counts and a panel holding both would compare the target
        against itself under another name (ratio ≡ 1.0, never trusting)."""
        seen, out = set(), []
        for m in (design.site.metal, *self.metals):
            if (el := element_symbol(m)) not in seen:
                seen.add(el)
                out.append(m)
        return out

    def verify(self, design: BinderDesign) -> Verdict:
        totals = metal_totals()
        target, *competitors = self.panel(design)
        if totals.get(element_symbol(target), 0) < self.min_metal_sites:
            return Verdict.defer(
                f"only {totals.get(element_symbol(target), 0)} precedented sites for {target} "
                f"(< {self.min_metal_sites}) — no basis to judge what its sites look like"
            )
        # a thin *competitor* is a gap in the reference table, not a fault in the design — drop it
        # from the panel rather than deferring a design the user never asked about it.
        thin = [m for m in competitors if totals.get(element_symbol(m), 0) < self.min_metal_sites]
        competitors = [m for m in competitors if m not in thin]
        if not competitors:
            return Verdict.defer(
                "no competitor metal with enough precedent to compare against"
                + (f" (dropped {', '.join(thin)})" if thin else "")
            )

        enr = motif_enrichment(design.site, tuple(competitors))
        if not enr.eligible:  # nobody is observed on this donor set often enough to own it
            return Verdict.defer(
                f"donor set {enr.donors} has < {MIN_MOTIF_HITS} precedents for every metal in the "
                "panel — too little evidence to say whose it is"
            )

        ratio, pct = enr.ratio, {m: f"{e * 100:.1f}%" for m, e in enr.enrichment.items()}
        reason = (
            f"{enr.donors} is most characteristic of {enr.preferred} "
            f"(target {enr.target} {pct[enr.target]}, ×{ratio:.2f} vs best competitor)"
            + (f"; dropped {', '.join(thin)} — too few precedents" if thin else "")
        )
        return Verdict(
            score=float(ratio / (1.0 + ratio)),  # 0.5 at parity, saturating
            trust=enr.preferred == enr.target and ratio >= self.trust_ratio,
            ood=False,
            reason=reason,
            metrics={
                "donors": enr.donors,
                "preferred": enr.preferred,
                "ratio": round(ratio, 3),
                "hits": enr.hits,
                "enrichment": {m: round(e, 5) for m, e in enr.enrichment.items()},
                "background": round(enr.background, 5),
            },
        )


class PrecedentVerifier:
    """Trusts a coordination motif with ≥ `min_hits` precedents; defers an unprecedented one
    (off the known manifold) or when the precedent search can't run. Defaults to the open
    MetalPDB search; pass a searcher for CSD-CrossMiner or a custom source."""

    def __init__(
        self,
        search: Callable[[CoordinationSite], PrecedentHits] | None = None,
        *,
        min_hits: int = 5,
    ):
        self.search = search or metalpdb_precedent_search
        self.min_hits = min_hits

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            hits = self.search(design.site)
        except Exception as e:  # search / data / licence failure ⇒ can't judge
            return Verdict.defer(f"precedent search failed: {type(e).__name__}")

        n = hits.nhits
        reason = f"{n} precedent(s) for {hits.motif}"
        if n == 0:  # unprecedented coordination — off the known manifold
            return Verdict.defer(reason + " — unprecedented")
        score = n / (n + self.min_hits)  # saturating: 0.5 at min_hits
        return Verdict(score, trust=n >= self.min_hits, ood=False, reason=reason)
