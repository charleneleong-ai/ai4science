"""CSD precedent check via CSD-CrossMiner (CCDC) — has this coordination motif been
seen before?

A cheap, strong prior over the CSD's >1.3M structures: a metal-donor environment with
many CSD matches is well-precedented (trust); an unprecedented one (no hits) sits off
the known manifold and is deferred for deeper physics scrutiny. As a *cheap* gate (a
CSD lookup, no GPU) it belongs at the front of `cascade`, so the expensive tiers only
run on designs with structural precedent.

The searcher is pluggable (inject one for tests); the default lazily drives CCDC's
CSD-CrossMiner, so it is licence-gated (CSD Python API) — its query construction should
be confirmed against the installed API version before quantitative use.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from ..core import BinderDesign, CoordinationSite, Verdict, element_symbol


@dataclass
class PrecedentHits:
    """CSD-CrossMiner support for one coordination motif."""

    nhits: int  # CSD structures matching the metal–donor environment
    motif: str  # e.g. "Ni-N3O2"


def _motif(site: CoordinationSite) -> str:
    """A compact label for a coordination environment, e.g. 'Ni-N3O2'."""
    from collections import Counter

    donors = "".join(f"{el}{n}" for el, n in sorted(Counter(site.ligand_elems).items()))
    return f"{element_symbol(site.metal)}-{donors}"


def _crossminer_search(site: CoordinationSite) -> PrecedentHits:
    """Search CSD-CrossMiner for the coordination motif. Requires `ccdc` + a licence;
    the pharmacophore/connectivity query should be confirmed against the installed API."""
    from ccdc.search import SubstructureSearch  # noqa: F401 — licence probe

    raise NotImplementedError(
        f"wire the CrossMiner query for {_motif(site)} against the licensed CSD API"
    )


class PrecedentVerifier:
    """Trusts a coordination motif with ≥ `min_hits` CSD precedents; defers an
    unprecedented one (off the known manifold) or when CrossMiner can't run."""

    def __init__(
        self,
        search: Callable[[CoordinationSite], PrecedentHits] | None = None,
        *,
        min_hits: int = 5,
    ):
        self._search = search or _crossminer_search
        self.min_hits = min_hits

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            hits = self._search(design.site)
        except Exception as e:  # search / licence failure ⇒ can't judge
            return Verdict.defer(f"CrossMiner search failed: {type(e).__name__}")

        n = hits.nhits
        reason = f"{n} CSD precedent(s) for {hits.motif}"
        if n == 0:  # unprecedented coordination — off the known manifold
            return Verdict.defer(reason + " — unprecedented")
        score = n / (n + self.min_hits)  # saturating: 0.5 at min_hits
        return Verdict(score, trust=n >= self.min_hits, ood=False, reason=reason)
