import numpy as np
from mindlessgen.symmetrization.rotation import CnRotation  # type: ignore
from mindlessgen.molecules.molecule import Molecule, ati_to_atlist  # type: ignore
from mindlessgen.prog.config import SymmetrizationConfig  # type: ignore


def test_rotation_modify_structure():
    # Create a mock molecule
    mol = Molecule()
    mol.xyz = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = ati_to_atlist(mol.ati)

    # Create a mock config
    config = SymmetrizationConfig()
    config.rotation = 2

    # Initialize the rotation class
    rotation = CnRotation(config)

    # Modify the structure
    modified_molecule = rotation.get_symmetric_structure(mol)

    # Check the modified molecule
    expected_xyz = np.array(
        [
            [2.5, 0.0, 0.0],
            [1.5, 1.0, 0.0],
            [-2.5, 0.0, 0.0],
            [-1.5, -1.0, 0.0],
        ]
    )
    # assert arrays are almost equal
    assert np.allclose(modified_molecule.xyz, expected_xyz, atol=1e-10)
    assert modified_molecule.num_atoms == 4
    assert np.array_equal(modified_molecule.ati, np.array([1, 1, 1, 1]))
    assert modified_molecule.charge == 0
    assert modified_molecule.uhf == 0
    ref_modified_atlist = np.zeros(103, dtype=int)
    ref_modified_atlist[1] = 4
    assert np.array_equal(modified_molecule.atlist, ref_modified_atlist)
