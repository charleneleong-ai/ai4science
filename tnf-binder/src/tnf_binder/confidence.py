"""Reading Boltz confidence output with the binder/target split made explicit.

The implementation now lives in `touchstone.boltz`, upstreamed because nothing about it is
TNF-specific: it reads chain-pair confidence out of a Boltz JSON and knows nothing about what
the chains are. Keeping a second copy here would repeat the mistake the composition verifier
already made, where the two drifted apart on a threshold boundary.

The two rules it encodes, both of which cost this campaign a wrong ranking before they were
written down: global ipTM is unusable on a multi-chain target because the target's own chain
pairs dominate it (+0.405 overstatement measured here), and diffusion samples aggregate on the
mean, because best-of-n ranks on sampling variance.
"""

from __future__ import annotations

from touchstone.boltz import (
    BINDER_CHAIN,
    BoltzConfidence,
    DesignConfidence,
    mean_or_zero,
    rank_designs,
)

__all__ = [
    "BINDER_CHAIN",
    "BoltzConfidence",
    "DesignConfidence",
    "mean_or_zero",
    "rank_designs",
]
