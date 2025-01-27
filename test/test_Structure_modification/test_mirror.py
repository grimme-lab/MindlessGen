import pytest
import numpy as np
from mindlessgen.Structure_modification.mirror import Mirror  # type: ignore
from mindlessgen.molecules.molecule import Molecule  # type: ignore
from mindlessgen.prog.config import StructureModConfig  # type: ignore


@pytest.fixture
def molecule():
    mol = Molecule()
    mol.xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.ndarray(103, dtype=int)
    mol.atlist[0] = 2
    return mol


@pytest.fixture
def config():
    return StructureModConfig()


def test_mirror_modify_structure(molecule, config):
    mirror = Mirror()
    modified_molecule = mirror.modify_structure(molecule, config)

    assert modified_molecule.num_atoms == 4
    assert modified_molecule.charge == molecule.charge * 2
    assert modified_molecule.uhf == molecule.uhf * 2
    assert np.array_equal(
        modified_molecule.ati, np.hstack((molecule.ati, molecule.ati))
    )
    assert np.array_equal(modified_molecule.xyz[:2], molecule.xyz)
    assert np.array_equal(
        modified_molecule.xyz[2:], np.array([[-0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    )
