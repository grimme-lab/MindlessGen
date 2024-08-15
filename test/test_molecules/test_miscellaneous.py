"""
Test the miscellaneous functions in the molecules module.
"""

import pytest
import numpy as np
from mindlessgen.molecules.miscellaneous import set_random_charge  # type: ignore


# CAUTION: We use 0-based indexing for atoms and molecules!
@pytest.mark.parametrize(
    "atom_types, expected_charge",
    [
        (np.array([5, 7, 0]), [-1, 1]),
        (np.array([5, 7, 0, 0]), [-2, 0, 2]),
    ],
    ids=["odd", "even"],
)
def test_set_random_charge(atom_types, expected_charge):
    """
    Test the set_random_charge function.
    """
    assert set_random_charge(atom_types) in expected_charge
