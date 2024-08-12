"""
Test the squaring function.
"""

from __future__ import annotations

import numpy as np

from mlmgen.molecules import generate_molecule


# @pytest.mark.parametrize("value", [1.0, 2, -3.0])
# def test_squarer(value: int | float) -> None:
#     expected = value * value
#     actual = generator(value)
#
#     assert pytest.approx(expected) == actual


def test_generate_molecule() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    natoms = generate_molecule()

    # assert that the array has this shape: natoms = np.zeros(87, dtype=int)
    assert natoms.shape == (87,)
    assert natoms.dtype == int
    assert np.sum(natoms) > 0


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
