import numpy as np
from mindlessgen.Structure_modification.inversion import Inversion  # type: ignore
from mindlessgen.molecules.molecule import Molecule  # type: ignore
from mindlessgen.prog.config import StructureModConfig  # type: ignore


def test_inversion_modify_structure():
    # Create a mock molecule
    mol = Molecule()
    mol.xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 2])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.ndarray(103, dtype=int)
    mol.atlist[0] = 1
    mol.atlist[7] = 1

    # Create a mock config
    config = StructureModConfig()

    # Initialize the Inversion class
    inversion = Inversion()

    # Perform the inversion
    modified_molecule = inversion.modify_structure(mol, config)

    # Check the modified molecule
    expected_xyz = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0]]
    )
    assert np.array_equal(modified_molecule.xyz, expected_xyz)
    assert modified_molecule.num_atoms == 4
    assert np.array_equal(modified_molecule.ati, np.array([1, 2, 1, 2]))
    assert modified_molecule.charge == 0
    assert modified_molecule.uhf == 0
