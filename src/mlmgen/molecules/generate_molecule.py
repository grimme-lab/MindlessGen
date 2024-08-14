"""
This module generates a random molecule with a random number of atoms.
"""

import copy
import numpy as np
from .molecule import Molecule
from .miscellaneous import set_random_charge


def generate_random_molecule(verbosity: int = 1) -> Molecule:
    """
    Generate a random molecule of type Molecule.
    """

    mol = Molecule()
    mol.atlist = generate_atom_list(verbosity)
    mol.num_atoms = np.sum(mol.atlist)
    mol.xyz, mol.ati = generate_coordinates(
        at=mol.atlist, scaling=3.0, dist_threshold=1.2, verbosity=verbosity
    )
    mol.charge = set_random_charge(mol.ati, verbosity)

    # if verbosity > 1, print the molecule
    if verbosity > 1:
        print(mol)

    return mol


def generate_atom_list(verbosity: int = 1) -> np.ndarray:
    """
    Generate a random molecule with a random number of atoms.
    """
    not_included = ()
    # not_included = (
    #     1,
    #     2,
    #     5,
    #     6,
    #     7,
    #     8,
    #     9,
    #     10,
    #     18,
    #     21,
    #     22,
    #     23,
    #     24,
    #     25,
    #     26,
    #     27,
    #     28,
    #     29,
    #     30,
    #     36,
    #     39,
    #     40,
    #     41,
    #     42,
    #     43,
    #     44,
    #     45,
    #     46,
    #     47,
    #     48,
    #     54,
    #     57,
    #     58,
    #     59,
    #     60,
    #     61,
    #     62,
    #     63,
    #     64,
    #     65,
    #     66,
    #     67,
    #     68,
    #     69,
    #     70,
    #     71,
    #     72,
    #     73,
    #     74,
    #     75,
    #     76,
    #     77,
    #     78,
    #     79,
    #     80,
    #     86,
    # )

    natoms = np.zeros(102, dtype=int)

    # Generating random atoms from whole PSE if no input file is found
    # Add random elements from the whole PSE
    # Define the number of atom types to be added
    numatoms_all = np.random.randint(1, 7)
    for _ in range(numatoms_all):
        # Define the atom type to be added
        ati = np.random.randint(0, 86)
        if verbosity > 1:
            print(f"Adding atom type {ati}...")
        while ati + 1 in not_included:
            ati = np.random.randint(0, 86)
            if verbosity > 1:
                print(f"Adding atom type {ati}...")
        # Add a random number of atoms of the defined type
        natoms[ati] = natoms[ati] + np.random.randint(0, 3)

    # > If too many alkaline and alkine earth metals are included, restart generation
    metals = (2, 3, 10, 11, 18, 19, 36, 37, 54, 55)
    nmetals = 0
    for i in metals:
        nmetals += natoms[i]
    if nmetals > 3:
        # reduce number of metals starting from 2, going to 55
        while nmetals > 3:
            for i in metals:
                if natoms[i] > 0:
                    natoms[i] = natoms[i] - 1
                    nmetals -= 1
                if nmetals <= 3:
                    break

    # If too many transition or lanthanide metals are included, reduce their number
    for i in get_metal_z():
        if natoms[i] > 1:
            natoms[i] = 1

    # Add Elements between B and F (5-9)
    for _ in range(5):
        i = np.random.randint(4, 10)
        natoms[i] = natoms[i] + np.random.randint(0, 3)

    # If no H is included, add H atoms
    if natoms[0] == 0:
        nat = np.sum(natoms)
        minnat = min(nat, 10)
        randint = np.random.rand()
        j = 1 + int(randint * minnat * 1.2)
        natoms[0] = natoms[0] + j

    return natoms


def generate_coordinates(
    at: np.ndarray, scaling: float, dist_threshold: float, verbosity: int = 1
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate random coordinates for a molecule.
    """

    EXPANSION_FACTOR = 1.3

    # eff_scaling is a deep copy of scaling
    eff_scaling = copy.deepcopy(scaling)
    xyz, ati = generate_random_coordinates(at)
    xyz = xyz * eff_scaling
    # do while check_distances is False
    while not check_distances(xyz, dist_threshold):
        if verbosity > 1:
            print(
                f"Distance check failed. Increasing expansion factor by {EXPANSION_FACTOR}..."
            )
        xyz, ati = generate_random_coordinates(at)
        eff_scaling = eff_scaling * EXPANSION_FACTOR
        xyz = xyz * eff_scaling

    return xyz, ati


def generate_random_coordinates(at: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate random coordinates for a molecule.
    """
    atilist: list[int] = []
    xyz = np.zeros((sum(at), 3))
    numatoms = 0
    for elem, count in enumerate(at):
        for m in range(count):
            # different rules for hydrogen
            if elem == 0:
                xyz[numatoms + m, :] = np.random.rand(3) * 3 - 1.5
            else:
                xyz[numatoms + m, :] = np.random.rand(3) * 2 - 1
            atilist.append(elem)

        numatoms += count

    ati = np.array(atilist, dtype=int)

    return xyz, ati


def check_distances(xyz: np.ndarray, threshold: float) -> bool:
    """
    Check if the distances between atoms are larger than a threshold.
    """
    # go through the atoms dimension of the xyz array
    for i in range(xyz.shape[0] - 1):
        for j in range(i + 1, xyz.shape[0]):
            r = np.linalg.norm(xyz[i, :] - xyz[j, :])
            if r < threshold:
                return False
    return True


def get_metal_z() -> list[int]:
    """
    Get the atomic numbers of transition metals and lanthanides, for which different rules apply.
    """
    metals = list(range(20, 30)) + list(range(38, 48)) + list(range(56, 80))

    return metals
