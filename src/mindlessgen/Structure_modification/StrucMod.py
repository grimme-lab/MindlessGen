"""
This module defines the abstract base class for all structure modification methods.
"""

from abc import ABC, abstractmethod

from mindlessgen.molecules.molecule import Molecule
from mindlessgen.prog.config import StructureModConfig


class StrucMod(ABC):
    """
    This abstract base class defines the interface for all structure modification methods.
    """

    def __init__(self, config: StructureModConfig):
        """
        Initialize the structure modification class.
        """
        self.cfg = config

    @abstractmethod
    def modify_structure(
        self,
        mol: Molecule,
    ) -> Molecule:
        """
        Define the structure modification process.

        Arguments:
        molecule (Molecule): Molecule to modify
        config (StructureModConfig): Configuration for the structure modification

        Returns:
        Molecule: Modified molecule
        """

    # All structure modification methods should be able to translate the molecule
    def translation(
        self,
        mol: Molecule,
    ) -> Molecule:
        """
        Translate the molecule and add a little extra distance to the x axis.
        """
        xyz = mol.xyz
        # Find the translation vector in x direction
        translation_vector = xyz.min(axis=0)[0]
        # missing distance to reach minimum x value
        xshift = self.cfg.distance / 2 - translation_vector
        for i in range(mol.num_atoms):
            xyz[i, 0] += xshift
        mol.xyz = xyz
        return mol
