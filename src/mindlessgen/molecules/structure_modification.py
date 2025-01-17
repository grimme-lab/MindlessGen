"""
Modifies the generated molecules to get symmetric NCI structures.
"""

import numpy as np
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
        xyz = mol.xyz
        translation_vector = xyz.min()
        print(f"Translation vector: {translation_vector}")
        for i in range(mol.num_atoms):
            if translation_vector < 0:
                xyz[i, 0] += -translation_vector + 0.5
        print(f"New xyz after translation: {xyz}")
        mol.xyz = xyz

    if config.mirror:
        print("Mirror the molecule")

        xyz_mirror = xyz.copy()
        for i in range(mol.num_atoms):
            xyz_mirror[i, 0] = -xyz[i, 0]

        xyz_combined = np.vstack((xyz, xyz_mirror))
        mol.xyz = xyz_combined
        mol.num_atoms += mol.num_atoms
        mol.ati = np.hstack((mol.ati, mol.ati))
        mol.charge += mol.charge
        mol.uhf += mol.uhf
        mol.atlist += mol.atlist
        print(f"New xyz after mirroring: {mol.xyz}")

        modified_molecule = mol.copy()
        print(modified_molecule.get_xyz_str())
        return modified_molecule
    if config.rotation:
        print("Rotating the molecule")
        n = config.n
        rotation_matrix = np.array(
            [
                [np.cos((2 * np.pi) / n), -np.sin((2 * np.pi) / n), 0],
                [np.sin((2 * np.pi) / n), np.cos((2 * np.pi) / n), 0],
                [0, 0, 1],
            ]
        )
        xyz_rotation = xyz.copy()
        for i in range(mol.num_atoms):
            xyz_rotation[i] = np.dot(rotation_matrix, xyz[i])

        print(f"New xyz after rotation: {xyz_rotation}")
        xyz_combined = np.vstack((xyz, xyz_rotation))
        mol.xyz = xyz_combined
        mol.num_atoms += mol.num_atoms
        mol.ati = np.hstack((mol.ati, mol.ati))
        mol.charge += mol.charge
        mol.uhf += mol.uhf
        mol.atlist += mol.atlist
        print(f"New xyz after mirroring: {mol.xyz}")

        modified_molecule = mol.copy()
        print(modified_molecule.get_xyz_str())
        return modified_molecule
    if config.inversion:
        print("Inverting the molecule")

        # inversion_matrix = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        xyz_inversion = xyz.copy()
        for i in range(mol.num_atoms):
            xyz_inversion[i, 0] = -xyz[i, 0]
            xyz_inversion[i, 1] = -xyz[i, 1]
            xyz_inversion[i, 2] = -xyz[i, 2]
        print(f"New xyz after inversion: {xyz_inversion}")

        xyz_combined = np.vstack((xyz, xyz_inversion))
        mol.xyz = xyz_combined
        mol.num_atoms += mol.num_atoms
        mol.ati = np.hstack((mol.ati, mol.ati))
        mol.charge += mol.charge
        mol.uhf += mol.uhf
        mol.atlist += mol.atlist
        print(f"New xyz after mirroring: {mol.xyz}")

        modified_molecule = mol.copy()
        print(modified_molecule.get_xyz_str())
        return modified_molecule

    print("structure_modification_mol is now working as a function")

    return mol
