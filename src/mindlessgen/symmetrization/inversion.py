"""
This class handles the inversion structure modification.
"""

import numpy as np

from .base import Symmetrizer
from ..molecules.molecule import Molecule


class Inversion(Symmetrizer):
    """
    This class handles the Inversion structure modification.
    """

    def _modify_structure(
        self,
        mol: Molecule,
        translation_distance: float,
    ) -> Molecule:
        """
        Invert the molecule.
        """
        mol = self.translation(mol, translation_distance)
        inversion_matrix = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        xyz_inversion = mol.xyz.copy()
        modified_molecule = mol.copy()
        for i in range(mol.num_atoms):
            xyz_inversion[i] = np.dot(inversion_matrix, mol.xyz[i])
        modified_molecule = self.combine_structure(
            mol, modified_molecule, xyz_inversion
        )
        return modified_molecule
