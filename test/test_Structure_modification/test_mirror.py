import pytest
import numpy as np
from mindlessgen.symmetrization.mirror import Mirror  # type: ignore
from mindlessgen.molecules.molecule import Molecule, ati_to_atlist  # type: ignore
from mindlessgen.prog.config import SymmetrizationConfig  # type: ignore


@pytest.fixture
def molecule():
    mol = Molecule()
    mol.xyz = np.array([[0.0, 4.0, 8.0], [2.0, 4.0, 8.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = ati_to_atlist(mol.ati)
    return mol


def test_mirror_modify_structure(molecule):
    # Create a mock config
    config = SymmetrizationConfig()
    config.distance = 4.0
    mirror = Mirror(config)
    modified_molecule = mirror.get_symmetric_structure(molecule)
    print(modified_molecule.xyz[:2])
    print(modified_molecule.xyz[2:])

    assert modified_molecule.num_atoms == 4
    assert modified_molecule.charge == molecule.charge * 2
    assert modified_molecule.uhf == molecule.uhf * 2
    assert np.array_equal(
        modified_molecule.ati, np.hstack((molecule.ati, molecule.ati))
    )
    assert np.array_equal(
        modified_molecule.xyz[:2], np.array([[2.0, 4.0, 8.0], [4.0, 4.0, 8.0]])
    )
    assert np.array_equal(
        modified_molecule.xyz[2:], np.array([[-2.0, 4.0, 8.0], [-4.0, 4.0, 8.0]])
    )
