"""
Test the miscellaneous functions in the molecules module.
"""

import pytest
import numpy as np
from mindlessgen.molecules.miscellaneous import set_random_charge  # type: ignore


# CAUTION: We use 0-based indexing for atoms and molecules!
@pytest.mark.parametrize(
    "atom_types, expected_charges, expected_uhf",
    [
        # Standard mode tests (no lanthanides)
        (np.array([4, 6, 0]), [-1, 1], 0),  # B, N, H (odd number of electrons)
        (
            np.array([4, 6, 0, 0]),
            [-2, 0, 2],
            0,
        ),  # B, N, H, H (even number of electrons)
        (np.array([5, 7, 0]), [-1, 1], 0),  # C, O, H (odd number of electrons)
        (
            np.array([5, 7, 0, 0]),
            [-2, 0, 2],
            0,
        ),  # C, O, H, H (even number of electrons)
        # Lanthanide mode tests
        (np.array([56, 7, 0]), [0], 0),  # La (57), O, H -> Lanthanide case, UHF = 0
        (
            np.array([59, 7, 7, 0]),
            [0],
            3,
        ),  # Nd (60), O, O, H -> Lanthanide case, UHF = 3
        (
            np.array([63, 6, 6]),
            [1],
            7,
        ),  # Gd (64), C, C -> Lanthanide case, UHF = 7, CHRG = 1
        (np.array([57, 7, 0]), [0], 1),  # Ce (58), O, H -> Lanthanide case, UHF = 1
        (np.array([69, 7, 0]), [0], 1),  # Yb (70), O, H -> Lanthanide case, UHF = 1
        (np.array([5, 7, 0, 98]), [0], 4),  # Es(99), O, C, H -> Actinides case, UHF = 4
        (
            np.array([59, 92, 0, 0]),
            [0],
            7,
        ),  # Nd(60), Np(93), H, H -> Lanthanide and Actinide case, UHF = 7
        (
            np.array([59, 91, 92, 5, 7, 0]),
            [0],
            10,
        ),  # Nd(60), U(92), Np(93), H, H -> Lanthanide and Actinides case, UHF = 10
    ],
    ids=[
        "B-N-H (standard, odd)",
        "B-N-H-H (standard, even)",
        "C-O-H (standard, odd)",
        "C-O-H-H (standard, even)",
        "La-O-H (lanthanide)",
        "Nd-O-O-H (lanthanide)",
        "Gd-C-N (lanthanide)",
        "Ce-O-H (lanthanide)",
        "Yb-O-H (lanthanide)",
        "Es-O-C-H (actinides)",
        "Nd-Np-H-H (lanthanide and actinide)",
        "Nd-Eu-Np-H-H (lanthanide and actinides)",
    ],
)
def test_set_random_charge(atom_types, expected_charges, expected_uhf):
    """
    Test the set_random_charge function for both standard and lanthanide modes.
    """
    charge, uhf = set_random_charge(atom_types, verbosity=1)
    assert charge in expected_charges
    assert uhf == expected_uhf
