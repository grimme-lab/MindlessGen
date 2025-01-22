"""
Test the structure modification class.
"""

import numpy as np

from mindlessgen.molecules.molecule import Molecule
from mindlessgen.molecules.structure_modification import structure_modification_mol
from mindlessgen.prog import StructureModConfig as config


def test_structure_modification_translation():
    config.mirror = False
    config.rotation = False
    config.inversion = False
    config.distance = 1.0
    mol = Molecule()
    mol.xyz = np.array([[-0.5, 0.0, 0.0], [0.5, 1.0, 1.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.zeros(103, dtype=int)
    mol.atlist[0] = 1
    modified_mol = structure_modification_mol(mol, config, verbosity=0)
    assert np.allclose(modified_mol.xyz, np.array([[0.5, 0.0, 0.0], [1.5, 1.0, 1.0]]))


def test_structure_modification_mirror():
    config.mirror = True
    config.rotation = False
    config.inversion = False
    mol = Molecule()
    mol.xyz = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.zeros(103, dtype=int)
    mol.atlist[0] = 1
    modified_mol = structure_modification_mol(mol, config, verbosity=0)
    expected_xyz = np.array(
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [-1.0, 1.0, 1.0]]
    )
    assert np.allclose(modified_mol.xyz, expected_xyz)


def test_structure_modification_rotation():
    config.mirror = False
    config.rotation = True
    config.n_fold_rotation = 4
    config.inversion = False
    mol = Molecule()
    mol.xyz = np.array([[1.0, 0.0, 0.0]])
    mol.num_atoms = 1
    mol.ati = np.array([0])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.zeros(103, dtype=int)
    mol.atlist[0] = 1
    modified_mol = structure_modification_mol(mol, config, verbosity=0)
    expected_xyz = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
    )
    assert np.allclose(modified_mol.xyz, expected_xyz)


def test_structure_modification_inversion():
    config.mirror = False
    config.rotation = False
    config.inversion = True
    mol = Molecule()
    mol.xyz = np.array([[1.0, 1.0, 1.0]])
    mol.num_atoms = 1
    mol.ati = np.array([1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.zeros(103, dtype=int)
    mol.atlist[0] = 1
    modified_mol = structure_modification_mol(mol, config, verbosity=0)
    expected_xyz = np.array([[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])
    assert np.allclose(modified_mol.xyz, expected_xyz)
