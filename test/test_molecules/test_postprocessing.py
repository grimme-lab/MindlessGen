"""
Test the postprocessing functions in the molecules module.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest
from mlmgen.molecules import detect_fragments  # type: ignore
from mlmgen.molecules.molecule import Molecule  # type: ignore


@pytest.fixture
def mol_H6O2B2Ne2I1Os1Tl1() -> Molecule:
    """
    Load the molecule H6O2B2Ne2I1Os1Tl1 from 'fixtures/H6O2B2Ne2I1Os1Tl1.xyz'.
    """
    mol = Molecule("H6O2B2Ne2I1Os1Tl1")
    # get the Path of this file
    path = Path(__file__).resolve().parent
    xyz_file = path / "fixtures/H6O2B2Ne2I1Os1Tl1.xyz"
    mol.read_xyz_from_file(xyz_file)
    return mol


@pytest.fixture
def mol_H2O2B2I1Os1() -> Molecule:
    """
    Load the molecule H2O2B2I1Os1 from 'fixtures/H2O2B2I1Os1.xyz'.
    """
    mol = Molecule("H2O2B2I1Os1")
    # get the Path of this file
    path = Path(__file__).resolve().parent
    xyz_file = path / "fixtures/H2O2B2I1Os1.xyz"
    mol.read_xyz_from_file(xyz_file)
    return mol


def test_detect_fragments_H6O2B2Ne2I1Os1Tl1(
    mol_H6O2B2Ne2I1Os1Tl1: Molecule,
    mol_H2O2B2I1Os1: Molecule,
) -> None:
    """
    Test the detection of fragments in the molecule H2O2B2I1Os1.
    """
    fragmols = detect_fragments(mol_H6O2B2Ne2I1Os1Tl1, verbosity=0)
    assert len(fragmols) == 7

    # check that the first fragment is equal to the molecule H2O2B2I1Os1
    assert np.allclose(fragmols[0].xyz, mol_H2O2B2I1Os1.xyz)
