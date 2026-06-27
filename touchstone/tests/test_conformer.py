import numpy as np
import pytest

from touchstone import conformer_seeds
from touchstone.core import BinderDesign, CoordinationSite


def _design() -> BinderDesign:
    site = CoordinationSite("Ni2+", np.zeros(3), np.empty((0, 3)), ())
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5)


def test_pluggable_generator_seeds_n():
    seeds = conformer_seeds(_design(), n=3, generate=lambda _d, k: ["conf"] * k)
    assert seeds == ["conf", "conf", "conf"]  # n forwarded to the generator


def test_default_is_licence_gated():
    with pytest.raises((ImportError, NotImplementedError)):
        conformer_seeds(_design())  # the default drives ccdc (not installed here)
