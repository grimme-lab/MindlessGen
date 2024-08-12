"""
Test the squaring function.
"""

from __future__ import annotations

import numpy as np
import pytest

from mlmgen.molecules import generate_molecule  # type: ignore
from mlmgen.molecules.molecule import Molecule  # type: ignore


def test_generate_molecule() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    mol = generate_molecule()

    assert mol.num_atoms > 0
    assert mol.num_atoms == np.sum(mol.atlist)
    assert mol.num_atoms == len(mol.xyz)
    # assert that sum of absolute value of mol.xyz is greater than 0
    assert np.sum(np.abs(mol.xyz)) > 0


@pytest.mark.parametrize(
    "value, expected_exception",
    [
        (3, None),  # Valid value
        (-1, ValueError),  # Invalid value
        ("three", TypeError),  # Invalid type
    ],
)
def test_num_atoms_property(value, expected_exception):
    mol = Molecule()
    if expected_exception:
        with pytest.raises(expected_exception):
            mol.num_atoms = value
    else:
        mol.num_atoms = value
        assert mol.num_atoms == value


@pytest.mark.parametrize(
    "value, expected_exception",
    [
        (1, None),  # Valid value
        ("positive", TypeError),  # Invalid type
    ],
)
def test_charge_property(value, expected_exception):
    mol = Molecule()
    if expected_exception:
        with pytest.raises(expected_exception):
            mol.charge = value
    else:
        mol.charge = value
        assert mol.charge == value


@pytest.mark.parametrize(
    "value, expected_exception",
    [
        (2, None),  # Valid value
        ("two", TypeError),  # Invalid type
    ],
)
def test_uhf_property(value, expected_exception):
    mol = Molecule()
    if expected_exception:
        with pytest.raises(expected_exception):
            mol.uhf = value
    else:
        mol.uhf = value
        assert mol.uhf == value


@pytest.mark.parametrize(
    "xyz_value, expected_exception",
    [
        (np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), None),  # Valid array
        (
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            TypeError,
        ),  # Invalid type
        (
            np.array([[0.0, 0.0]]),
            ValueError,
        ),  # Invalid shape
    ],
)
def test_xyz_property(xyz_value, expected_exception):
    mol = Molecule()
    if expected_exception:
        with pytest.raises(expected_exception):
            mol.xyz = xyz_value
    else:
        mol.num_atoms = xyz_value.shape[0]  # Set num_atoms based on xyz array shape
        mol.xyz = xyz_value
        np.testing.assert_array_equal(mol.xyz, xyz_value)


def test_ati_property():
    mol = Molecule()
    # create an empty array with dimension 1 and length 86
    mol.atlist = np.zeros(102, dtype=int)
    assert mol.atlist.shape == (102,)

    # set the first element to 4 and the sixth element to 2
    mol.atlist[0] = 4
    mol.atlist[5] = 2

    # generate the sum formula, which should be 'C4H2'
    assert mol.sum_formula() == "C2H4"


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
