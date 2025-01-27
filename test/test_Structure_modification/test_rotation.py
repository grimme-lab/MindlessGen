import numpy as np
from mindlessgen.Structure_modification.rotation import CnRotation  # type: ignore
from mindlessgen.molecules.molecule import Molecule  # type: ignore
from mindlessgen.prog.config import StructureModConfig  # type: ignore


def test_modify_structure():
    # Create a mock molecule
    mol = Molecule()
    mol.xyz = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    mol.num_atoms = 2
    mol.ati = np.array([1, 1])
    mol.charge = 0
    mol.uhf = 0
    mol.atlist = np.ndarray(103, dtype=int)
    mol.atlist[0] = 2

    # Create a mock config
    config = StructureModConfig()
    config.rotation = 2

    # Initialize the rotation class
    rotation = CnRotation()

    # Modify the structure
    modified_molecule = rotation.modify_structure(mol, config)

    # Check the modified molecule
    assert modified_molecule.num_atoms == 4
    assert modified_molecule.ati.shape == (4,)
    assert modified_molecule.xyz.shape == (4, 3)
    assert modified_molecule.charge == 0
    assert modified_molecule.uhf == 0
    assert modified_molecule.atlist[0] == 4
