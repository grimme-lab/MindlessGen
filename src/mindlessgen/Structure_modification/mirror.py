import numpy as np

from .StrucMod import StrucMod, Translation
from ..prog.config import StructureModConfig
from ..molecules.molecule import Molecule


class Mirror(StrucMod):
    """
    This class handles the translation structure modification.
    """

    def __init__(self):
        """
        Initialize the symmetrie operation structure modification class.
        """
        self.cfg = StructureModConfig()

    def modify_structure(
        self,
        mol: Molecule,
        translation: Translation,
        config: StructureModConfig,
    ) -> Molecule:
        """
        Mirror the molecule.
        """

        mol = translation.translation(mol, config)
        xyz = mol.xyz
        mirror_matrix = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        xyz_mirror = xyz.copy()
        modified_molecule = mol.copy()
        for i in range(mol.num_atoms):
            xyz_mirror[i] = np.dot(mirror_matrix, xyz[i])
        # Combine the original and the mirror image
        xyz_combined = np.vstack((xyz, xyz_mirror))
        modified_molecule.xyz = xyz_combined
        modified_molecule.num_atoms += mol.num_atoms
        modified_molecule.ati = np.hstack((modified_molecule.ati, mol.ati))
        modified_molecule.charge += mol.charge
        modified_molecule.uhf += mol.uhf
        modified_molecule.atlist += mol.atlist
        modified_molecule.set_name_from_formula()
        return modified_molecule
