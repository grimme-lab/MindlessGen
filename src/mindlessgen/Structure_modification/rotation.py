"""
This class handles the rotation structure modification.
"""

import numpy as np

from .StrucMod import StrucMod
from ..molecules.molecule import Molecule


class CnRotation(StrucMod):
    """
    This class handles the rotation structure modification.
    """

    def modify_structure(
        self,
        mol: Molecule,
    ) -> Molecule:
        """
        Rotate the molecule around the z-axis.
        """
        mol = self.translation(mol)
        xyz = mol.xyz
        n = self.cfg.rotation
        rotation_matrix = np.array(
            [
                [np.cos((2 * np.pi) / n), -np.sin((2 * np.pi) / n), 0],
                [np.sin((2 * np.pi) / n), np.cos((2 * np.pi) / n), 0],
                [0, 0, 1],
            ]
        )
        xyz_rotation = xyz.copy()
        modified_molecule = mol.copy()
        # For loop to get n times the molecule
        for _ in range(1, n):
            for k in range(mol.num_atoms):
                xyz_rotation[k] = np.dot(rotation_matrix, xyz_rotation[k])
            xyz = np.vstack((xyz, xyz_rotation))
            modified_molecule.num_atoms += mol.num_atoms
            modified_molecule.ati = np.hstack((modified_molecule.ati, mol.ati))
            modified_molecule.charge += mol.charge
            modified_molecule.uhf += mol.uhf
            modified_molecule.atlist += mol.atlist
        modified_molecule.xyz = xyz
        modified_molecule.set_name_from_formula()
        return modified_molecule
