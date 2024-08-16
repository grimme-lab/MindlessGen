"""
Test the squaring function.
"""

from __future__ import annotations
import pytest
import numpy as np
from mindlessgen.molecules import (  # type: ignore
    generate_random_molecule,
    generate_coordinates,
    generate_atom_list,
    check_distances,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
    get_alkali_metals,
    get_alkaline_earth_metals,
)
from mindlessgen.molecules.molecule import Molecule  # type: ignore
from mindlessgen.prog import ConfigManager  # type: ignore


def test_generate_molecule() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    # create a ConfigManager object with verbosity set to 0
    config = ConfigManager()
    config.general.verbosity = 0
    mol = generate_random_molecule(config.generate, config.general.verbosity)

    assert mol.num_atoms > 0
    assert mol.num_atoms == np.sum(mol.atlist)
    assert mol.num_atoms == len(mol.xyz)
    # assert that sum of absolute value of mol.xyz is greater than 0
    assert np.sum(np.abs(mol.xyz)) > 0


def test_generate_coordinates() -> None:
    """
    Test the generation of coordinates.
    """
    mol = Molecule()
    # create an empty array with dimension 1 and length 86
    mol.atlist = np.zeros(102, dtype=int)
    assert mol.atlist.shape == (102,)

    # set the first element to 4 and the sixth element to 2
    mol.atlist[0] = 4
    mol.atlist[5] = 2

    mol.xyz, mol.ati = generate_coordinates(mol.atlist, 3.0, 1.2)
    # assert that the shape of mol.xyz is (6, 3)
    assert mol.xyz.shape == (6, 3)
    assert mol.num_atoms == 6
    assert mol.num_atoms == np.sum(mol.atlist)
    assert mol.num_atoms == len(mol.xyz)
    assert mol.num_atoms == len(mol.ati)


def test_generate_atom_list() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    atlist = generate_atom_list()
    assert atlist.shape == (102,)
    assert np.sum(atlist) > 0
    # check that for the transition and lanthanide metals, the occurence is never greater than 3
    all_metals = (
        get_three_d_metals()
        + get_four_d_metals()
        + get_five_d_metals()
        + get_lanthanides()
    )
    for z in all_metals:
        assert atlist[z] <= 3
    alkmetals = get_alkali_metals() + get_alkaline_earth_metals()
    # check that the sum of alkali and alkaline earth metals is never greater than 3
    assert np.sum([atlist[z] for z in alkmetals]) <= 3


@pytest.mark.parametrize(
    "xyz, threshold, expected, description",
    [
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            0.5,
            True,
            "Two atoms with distance greater than threshold (1.0 > 0.5)",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]]),
            0.5,
            False,
            "Two atoms with distance less than threshold (0.4 < 0.5)",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            0.5,
            True,
            "Three atoms in a line with distances greater than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            0.5,
            False,
            "Three atoms with one pair close together: distance between first two is less than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0]]),
            0.5,
            True,
            "Single atom, no distances to compare",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            0.5,
            False,
            "Two atoms at identical positions: distance is zero, less than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            1.7320,
            True,
            "Two atoms with diagonal distance just above threshold (sqrt(3) ≈ 1.732)",
        ),
    ],
    ids=[
        "far_apart",
        "close_together",
        "three_in_line",
        "three_with_one_close",
        "single_atom",
        "two_identical",
        "diagonal_distance",
    ],
)
def test_check_distances(xyz, threshold, expected, description):
    assert check_distances(xyz, threshold) == expected


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
