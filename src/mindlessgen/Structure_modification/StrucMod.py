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

    def __init__(self):
        """
        Initialize the structure modification class.
        """
        self.cfg = StructureModConfig()

    @abstractmethod
    def modify_structure(
        self,
        mol: Molecule,
        config: StructureModConfig,
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
        config: StructureModConfig,
    ) -> Molecule:
        """
        Translate the molecule and add a little extra distance to the x axis.
        """
        xyz = mol.xyz.copy()
        translation_vector = xyz.min()
        for i in range(mol.num_atoms):
            if translation_vector < 0:
                xyz[i, 0] += (
                    -translation_vector + config.distance / 2
                )  # extra distance to avoid that an atom is at 0
        mol.xyz = xyz
        return mol
