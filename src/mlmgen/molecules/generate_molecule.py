"""
This module generates a random molecule with a random number of atoms.
"""

import copy

import numpy as np

from .molecule import Molecule


def generate_molecule(verbosity: int = 1) -> Molecule:
    """
    Generate a random molecule of type Molecule.
    """

    mol = Molecule()
    mol.ati = generate_atom_list(verbosity)
    mol.num_atoms = np.sum(mol.ati)
    mol.xyz = generate_coordinates(mol.ati, 3.0, 1.2)

    # if verbosity > 0, print the molecule and its sum formula
    if verbosity > 0:
        print(mol)
        print(mol.sum_formula())

    return mol


def generate_atom_list(verbosity: int = 1) -> np.ndarray:
    """
    Generate a random molecule with a random number of atoms.
    """
    not_included = (
        1,
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        18,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        36,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        54,
        57,
        58,
        59,
        60,
        61,
        62,
        63,
        64,
        65,
        66,
        67,
        68,
        69,
        70,
        71,
        72,
        73,
        74,
        75,
        76,
        77,
        78,
        79,
        80,
        86,
    )

    natoms = np.zeros(102, dtype=int)

    # Generating random atoms from whole PSE if no input file is found
    numatoms_all = np.random.randint(0, 3)
    for _ in range(1, numatoms_all + 1):
        ati = np.random.randint(0, 85)
        while (ati + 1) in not_included:
            ati = np.random.randint(0, 85)
        natoms[ati] = natoms[ati] + np.random.randint(0, 3)

    # Add Elements between B and F (5-9)
    for _ in range(1, 5):
        i = np.random.randint(4, 9)
        natoms[i] = natoms[i] + np.random.randint(0, 3)

    # If no H is included, add H atoms
    if natoms[0] == 0:
        natoms[0] = np.random.randint(0, 3)

    # # > If too many metals are included, restart generation
    # metals = (3, 4, 11, 12, 19, 20, 37, 38, 55, 56)
    # nmetals = 0
    # for i in metals:
    #     if natoms[i] > 0:
    #         nmetals += 1
    # if nmetals > 3:
    #     print("Too many metals. Restarting...")
    #     generate_molecule(natoms)

    return natoms


def generate_coordinates(
    at: np.ndarray, scaling: float, dist_threshold: float
) -> np.ndarray:
    """
    Generate random coordinates for a molecule.
    """

    # eff_scaling is a deep copy of scaling
    eff_scaling = copy.deepcopy(scaling)
    xyz = generate_random_coordinates(at)
    xyz = xyz * eff_scaling
    # do while check_distances is False
    while not check_distances(xyz, dist_threshold):
        xyz = generate_random_coordinates(at)
        eff_scaling = eff_scaling * 1.3
        xyz = xyz * eff_scaling

    return xyz


def generate_random_coordinates(at: np.ndarray) -> np.ndarray:
    """
    Generate random coordinates for a molecule.
    """
    xyz = np.zeros((sum(at), 3))
    numatoms = 0
    for k, i in enumerate(at):
        for m in range(i):
            # different rules for hydrogen
            if k == 0:
                xyz[numatoms + m, :] = np.random.rand(3) * 3 - 1.5
            else:
                xyz[numatoms + m, :] = np.random.rand(3) * 2 - 1
        numatoms += i

    return xyz


def check_distances(xyz: np.ndarray, threshold: float) -> bool:
    """
    Check if the distances between atoms are larger than a threshold.
    """
    # go through the atoms dimension of the xyz array
    for i in range(1, xyz.shape[1]):
        for j in range(i - 1):
            r = np.linalg.norm(xyz[:, i] - xyz[:, j])
            if r < threshold:
                return False
    return True
