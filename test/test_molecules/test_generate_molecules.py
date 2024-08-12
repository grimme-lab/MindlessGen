"""
Test the squaring function.
"""

from __future__ import annotations

import numpy as np

from mlmgen.molecules import generate_molecule, generate_coordinates  # type: ignore
from mlmgen.molecules.molecule import Molecule  # type: ignore


def test_generate_molecule() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    mol = generate_molecule(verbosity=0)

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


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
