"""
Test the squaring function.
"""

from __future__ import annotations

import numpy as np
import pytest

from mlmgen.mymath import square_a_number


@pytest.mark.parametrize("value", [1.0, 2, -3.0])
def test_squarer(value: int | float) -> None:
    expected = value * value
    actual = square_a_number(value)

    assert pytest.approx(expected) == actual


def test_squarer_fail() -> None:
    with pytest.raises(TypeError):
        square_a_number("2")  # type: ignore


def test_dummy() -> None:
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
