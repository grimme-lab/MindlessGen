"""
This module generates a random molecule with a random number of atoms.
"""

import numpy as np


def generate_molecule(verbosity: int | None = 1) -> np.ndarray:
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

    natoms = np.zeros(87, dtype=int)

    # Generating random atoms from whole PSE if no input file is found
    numatoms_all = np.random.randint(1, 4)
    for _ in range(1, numatoms_all + 1):
        ati = np.random.randint(1, 87)
        while ati in not_included:
            ati = np.random.randint(1, 87)
        natoms[ati] = natoms[ati] + np.random.randint(1, 4)

    # Add Elements between B and F (5-9)
    for _ in range(1, 5):
        i = np.random.randint(5, 10)
        natoms[i] = natoms[i] + np.random.randint(1, 4)

    # If no H is included, add H atoms
    if natoms[1] == 0:
        natoms[1] = np.random.randint(1, 4)

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
