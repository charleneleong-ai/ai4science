"""Coordination-motif precedent check — has this metal environment been seen before?

A cheap, strong prior: a metal–donor motif with many precedents is well-precedented (trust);
an unprecedented one (no hits) sits off the known manifold and is deferred for deeper scrutiny.
As a *cheap* gate (a table lookup, no GPU) it belongs at the front of `cascade`, so the
expensive tiers only run on designs with structural precedent.

The default searcher is **open**: `metalpdb_precedent_search` reads a bundled motif→count table mined
from MetalPDB (scripts/build_metalpdb_precedents.py) — the licence-free, metalloprotein-specific
analog of a CSD search. A CSD-CrossMiner searcher (`_crossminer_search`, licence-gated) is
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

_PRECEDENTS_DATA = Path(__file__).parent.parent / "data" / "metalpdb_precedents.json"


@dataclass
class PrecedentHits:
    """Precedent support for one coordination motif."""

    nhits: int  # structures matching the metal–donor environment
    motif: str  # e.g. "Ni-N3O2"


def _motif(site: CoordinationSite) -> str:
    """A compact label for a coordination environment, e.g. 'Ni-N3O2'."""
    donors = "".join(f"{el}{n}" for el, n in sorted(Counter(site.ligand_elems).items()))
    return f"{element_symbol(site.metal)}-{donors}"


@lru_cache(maxsize=1)
def _precedents_table() -> dict[str, int]:
    """The bundled MetalPDB motif→count table (built by scripts/build_metalpdb_precedents.py)."""
    return json.loads(_PRECEDENTS_DATA.read_text())


def metalpdb_precedent_search(site: CoordinationSite) -> PrecedentHits:
    """Open precedent via the bundled MetalPDB motif→count table — licence-free, the
    metalloprotein analog of a CSD-CrossMiner search. Raises if the table isn't bundled."""
    motif = _motif(site)
    return PrecedentHits(nhits=int(_precedents_table().get(motif, 0)), motif=motif)


def _crossminer_search(site: CoordinationSite) -> PrecedentHits:
    """Search CSD-CrossMiner for the coordination motif. Requires `ccdc` + a licence;
    the pharmacophore/connectivity query should be confirmed against the installed API."""
    from ccdc.search import SubstructureSearch  # noqa: F401 — licence probe

    raise NotImplementedError(
        f"wire the CrossMiner query for {_motif(site)} against the licensed CSD API"
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
        self._search = search or metalpdb_precedent_search
        self.min_hits = min_hits

    def verify(self, design: BinderDesign) -> Verdict:
        try:
            hits = self._search(design.site)
        except Exception as e:  # search / data / licence failure ⇒ can't judge
            return Verdict.defer(f"precedent search failed: {type(e).__name__}")

        n = hits.nhits
        reason = f"{n} precedent(s) for {hits.motif}"
        if n == 0:  # unprecedented coordination — off the known manifold
            return Verdict.defer(reason + " — unprecedented")
        score = n / (n + self.min_hits)  # saturating: 0.5 at min_hits
        return Verdict(score, trust=n >= self.min_hits, ood=False, reason=reason)
