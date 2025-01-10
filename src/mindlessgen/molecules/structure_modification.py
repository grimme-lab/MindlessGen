"""
Modifies the generated molecules to get symmetric NCI structures.
"""

from .molecule import Molecule
from ..prog import StructureModConfig


def structure_modification_mol(
    mol: Molecule, config: StructureModConfig, verbosity: int = 1
) -> Molecule:
    """
    Modifies the generated molecules to get symmetric NCI structures.

    Arguments:
    mol (Molecule): Molecule to modify

    Returns:
    Molecule: Modified molecule
    """

    if verbosity > 2:
        print("Modifying molecule structure...")

    if config.translation:
        print("Translating the molecule")
    if config.mirroring:
        print("Mirroring the molecule")
    if config.rotation:
        print("Rotating the molecule")
    if config.inversion:
        print("Inverting the molecule")

    print("structure_modification_mol is now working as a function")

    return mol
