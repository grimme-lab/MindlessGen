"""
Test the squaring function.
"""

from __future__ import annotations
import numpy as np
import pytest
from mindlessgen.molecules.molecule import Molecule  # type: ignore


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


def test_atlist_property():
    mol = Molecule()
    # create an empty array with dimension 1 and length 86
    mol.atlist = np.zeros(102, dtype=int)
    assert mol.atlist.shape == (102,)

    # set the first element to 4 and the sixth element to 2
    mol.atlist[0] = 4
    mol.atlist[5] = 2

    # generate the sum formula, which should be 'C4H2'
    assert mol.sum_formula() == "C2H4"


@pytest.mark.parametrize(
    "ati_value, num_atoms_value, expected_exception",
    [
        (np.array([1, 6]), 2, None),  # Valid array
        ([1, 6], 2, TypeError),  # Invalid type (list instead of numpy array)
        (np.array([1, 6, 8]), 2, ValueError),  # Invalid shape (length mismatch)
    ],
)
def test_ati_property(ati_value, num_atoms_value, expected_exception):
    mol = Molecule()
    mol.num_atoms = num_atoms_value
    if expected_exception:
        with pytest.raises(expected_exception):
            mol.ati = ati_value
    else:
        mol.ati = ati_value
        np.testing.assert_array_equal(mol.ati, ati_value)


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
