"""
This module defines the abstract base class for all structure modification methods.
"""

from abc import ABC, abstractmethod

from mindlessgen.molecules.molecule import Molecule
from mindlessgen.prog.config import StructureModConfig


class Translation(ABC):
    """
    This class handles the translation structure modification.
    """

    def translation(
        self,
        mol: Molecule,
        config: StructureModConfig,
    ) -> Molecule:
        """
        Translate the molecule and add a little extra distance to the x axis.
        """
        xyz = mol.xyz
        translation_vector = xyz.min()
        for i in range(mol.num_atoms):
            if translation_vector < 0:
                xyz[i, 0] += (
                    -translation_vector + config.distance / 2
                )  # extra distance to avoid that an atom is at 0
        mol.xyz = xyz
        return mol


class StrucMod(ABC):
    """
    This abstract base class defines the interface for all structure modification methods.
    """

    @abstractmethod
    def __init__(self):
        """
        Initialize the structure modification class.
        """
        self.cfg = StructureModConfig()

    @abstractmethod
    def modify_structure(
        self,
        molecule: Molecule,
        translation: Translation,
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
