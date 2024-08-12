"""
Test the squaring function.
"""

from __future__ import annotations

import numpy as np

from mlmgen.molecules import generate_molecule  # type: ignore


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


def test_dummy() -> None:
    """
    Test the dummy function.
    """
    # show effect of `conftest.py` by setting printoptions
    print(np.array([1.0 / 3.0]))
