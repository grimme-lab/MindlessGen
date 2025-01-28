import numpy as np
from mindlessgen.symmetrization.inversion import Inversion  # type: ignore
from mindlessgen.molecules.molecule import Molecule, ati_to_atlist  # type: ignore
from mindlessgen.prog.config import SymmetrizationConfig  # type: ignore


def test_inversion_modify_structure():
    # Create a mock molecule
    mol = Molecule()
    mol.xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 2])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = ati_to_atlist(mol.ati)

    # Create a mock config
    config = SymmetrizationConfig()
    config.distance = 3.0

    # Initialize the Inversion class
    inversion = Inversion(config)

    # Perform the inversion
    modified_molecule = inversion.get_symmetric_structure(mol)

    # Check the modified molecule
    expected_xyz = np.array(
        [[1.5, 2.0, 3.0], [4.5, 5.0, 6.0], [-1.5, -2.0, -3.0], [-4.5, -5.0, -6.0]]
    )
    assert np.array_equal(modified_molecule.xyz, expected_xyz)
    assert modified_molecule.num_atoms == 4
    assert np.array_equal(modified_molecule.ati, np.array([1, 2, 1, 2]))
    assert modified_molecule.charge == 0
    assert modified_molecule.uhf == 0
