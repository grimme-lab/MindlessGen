"""
This class handles the inversion structure modification.
"""

import numpy as np

from .base import Symmetrizer
from ..molecules.molecule import Molecule


class Mirror(Symmetrizer):
    """
    This class handles the translation structure modification.
    """

    def _modify_structure(
        self,
        mol: Molecule,
        translation_distance: float,
    ) -> Molecule:
        """
        Mirror the molecule.
        """
        mol = self.translation(mol, translation_distance)
        mirror_matrix = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        xyz_mirror = mol.xyz.copy()
        modified_molecule = mol.copy()
        for i in range(mol.num_atoms):
            xyz_mirror[i] = np.dot(mirror_matrix, mol.xyz[i])
        # Combine the original and the mirror image
        modified_molecule = self.combine_structure(mol, modified_molecule, xyz_mirror)
        return modified_molecule
