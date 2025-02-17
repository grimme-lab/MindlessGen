"""
This module defines the abstract base class for all symmetrization methods.
"""

from abc import ABC, abstractmethod

import numpy as np

from ..prog.config import SymmetrizationConfig
from ..molecules import check_distances, Molecule


class Symmetrizer(ABC):
    """
    This abstract base class defines the interface for all symmetrization methods.
    """

    def __init__(self, config: SymmetrizationConfig):
        """
        Initialize the symmetrization class.
        """
        self.cfg = config

    @abstractmethod
    def _modify_structure(
        self,
        mol: Molecule,
        translation_distance: float,
    ) -> Molecule:
        """
        Define the symmetrization process.

        Arguments:
        molecule (Molecule): Molecule to modify
        translation_distance (float): Distance to translate the molecule

        Returns:
        Molecule: Modified molecule
        """

    # All structure modification methods should be able to translate the molecule
    def translation(
        self,
        mol: Molecule,
        x_vec_translation: float | None = None,
    ) -> Molecule:
        """
        Translate the molecule and add a little extra distance to the x axis.

        Arguments:
        mol (Molecule): Molecule to translate
        x_vec_translation (float): Translation distance in x direction

        Returns:
        Molecule: Translated molecule
        """
        xyz = mol.xyz
        # Find the translation vector in x direction
        translation_vector = xyz.min(axis=0)[0]
        # missing distance to reach minimum x value
        if x_vec_translation is None:
            x_vec_translation = self.cfg.distance
        xshift = x_vec_translation / 2 - translation_vector
        for i in range(mol.num_atoms):
            xyz[i, 0] += xshift
        mol.xyz = xyz
        return mol

    def combine_structure(
        self,
        original_mol: Molecule,
        new_mol: Molecule,
        addxyz: np.ndarray,
    ) -> Molecule:
        """
        Combine the original and the modified molecule.

        Arguments:
        original_mol (Molecule): Original molecule
        new_mol (Molecule): Modified molecule
        addxyz (np.ndarray): Coordinates to add to the molecule

        Returns:
        Molecule: Combined molecule
        """
        new_mol.xyz = np.vstack((new_mol.xyz, addxyz))
        # add dimension of addxyz to the number of atoms
        new_mol.num_atoms += original_mol.num_atoms
        new_mol.ati = np.hstack((new_mol.ati, original_mol.ati))
        new_mol.charge += original_mol.charge
        new_mol.uhf += original_mol.uhf
        new_mol.atlist += original_mol.atlist
        new_mol.set_name_from_formula()
        return new_mol

    def get_symmetric_structure(
        self,
        mol: Molecule,
    ) -> Molecule:
        """
        Get the symmetric structure of the molecule.

        Arguments:
        mol (Molecule): Molecule to symmetrize

        Returns:
        Molecule: Symmetric molecule
        """
        orig_mol = mol.copy()
        if not check_distances(mol.xyz, mol.ati, 0.6):
            raise RuntimeError("The molecule has too close contacts.")
        translation_distance = self.cfg.distance
        mol = self._modify_structure(orig_mol, translation_distance)
        # scale the reference covalent radii for the check of interatomic distances by 0.8
        # to catch only REAL close contacts
        counter = 0
        while not check_distances(mol.xyz, mol.ati, scale_minimal_distance=0.6):
            counter += 1
            if counter > 1000:
                raise RuntimeError(
                    "Could not find symmetric structure without clashes."
                )
            translation_distance += 0.1
            mol = self._modify_structure(orig_mol, translation_distance)
        return mol
