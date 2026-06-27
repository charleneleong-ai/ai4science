"""CSD conformer seeding (CCDC ConformerGenerator) — physically-likely starting
geometries to cut MD sampling.

Rather than letting the expensive MD tier explore from one arbitrary pose, generate the
binder scaffold's CSD-likely conformers and seed MD from them — fewer, better-placed
starts converge faster. A utility, not a verifier.

Pluggable generator; the default lazily drives CCDC's ConformerGenerator (licence-gated).
Honest caveat: ConformerGenerator models *organic torsions*, not the metal coordination
sphere — it seeds the binder framework, not the first-shell geometry itself.
"""

from __future__ import annotations

from typing import Callable

from ..core import BinderDesign


def _ccdc_conformers(design: BinderDesign, n: int) -> list:
    from ccdc.conformer import ConformerGenerator  # noqa: F401 — licence probe

    raise NotImplementedError(
        "build a ccdc Molecule from the design's scaffold and run ConformerGenerator"
    )


def conformer_seeds(
    design: BinderDesign,
    *,
    n: int = 10,
    generate: Callable[[BinderDesign, int], list] | None = None,
) -> list:
    """Up to `n` CSD-likely conformer geometries to seed the MD tier from. `generate`
    is pluggable; the default drives CCDC's ConformerGenerator (licence-gated)."""
    return (generate or _ccdc_conformers)(design, n)
