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

    xyz = mol.xyz
    translation_vector = xyz.min()
    for i in range(mol.num_atoms):
        if translation_vector < 0:
            xyz[i, 0] += (
                -translation_vector + config.distance / 2
            )  # extra distance to avoid that an atom is at 0
    mol.xyz = xyz

    if config.mirror:
        mirror_matrix = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        xyz_mirror = xyz.copy()
        modified_molecule = mol.copy()
        for i in range(mol.num_atoms):
            xyz_mirror[i] = np.dot(mirror_matrix, xyz[i])

        # Combine the original and the mirror image
        xyz_combined = np.vstack((xyz, xyz_mirror))
        modified_molecule = update_molecule_objekt(modified_molecule, mol, xyz_combined)
        modified_molecule.set_name_from_formula()
        return modified_molecule
    elif config.rotation:
        n = config.n_fold_rotation
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
        for i in range(1, n):
            for k in range(mol.num_atoms):
                xyz_rotation[k] = np.dot(rotation_matrix, xyz_rotation[k])
            xyz = np.vstack((xyz, xyz_rotation))
            modified_molecule = update_molecule_objekt(modified_molecule, mol, xyz)

        modified_molecule.set_name_from_formula()
        return modified_molecule
    elif config.inversion:
        inversion_matrix = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        xyz_inversion = xyz.copy()
        modified_molecule = mol.copy()
        for i in range(mol.num_atoms):
            xyz_inversion[i] = np.dot(inversion_matrix, xyz[i])

        xyz_combined = np.vstack((xyz, xyz_inversion))
        modified_molecule = update_molecule_objekt(modified_molecule, mol, xyz_combined)
        modified_molecule.set_name_from_formula()
        return modified_molecule

    print("structure_modification_mol is now working as a function")

    return mol


def update_molecule_objekt(
    modified_molecule: Molecule, mol: Molecule, xyz: np.ndarray
) -> Molecule:
    """
    Updates the molecule object to the new twice as large molecule objekt.
    """
    modified_molecule.xyz = xyz
    modified_molecule.num_atoms += mol.num_atoms
    modified_molecule.ati = np.hstack((modified_molecule.ati, mol.ati))
    modified_molecule.charge += mol.charge
    modified_molecule.uhf += mol.uhf
    modified_molecule.atlist += mol.atlist
    return modified_molecule
