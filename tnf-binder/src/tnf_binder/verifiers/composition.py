"""Reject degenerate sequences, independently of any model score.

The implementation now lives in `touchstone.composition`, upstreamed because nothing about it
is TNF-specific: it judges `design.sequence` and no more. Keeping a second copy here was a
liability rather than a convenience — the two had already drifted on whether a sequence sitting
exactly at `max_single_fraction` is rejected, and a scoring bug once had to be fixed in both.

It exists because interface-confidence models score poly-alanine designs confidently. One such
sequence (59% alanine) reached 0.556 binder-to-target ipTM and sat at the top of a 40-design
shortlist. A composition check is cheap, model-free, and catches it immediately.
"""

from __future__ import annotations

from touchstone.composition import CompositionVerifier

__all__ = ["CompositionVerifier"]
