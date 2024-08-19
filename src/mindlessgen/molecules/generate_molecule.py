"""
This module generates a random molecule with a random number of atoms.
"""

import copy
import numpy as np
from ..prog import GenerateConfig
from .molecule import Molecule
from .miscellaneous import set_random_charge


def generate_random_molecule(
    config_generate: GenerateConfig, verbosity: int
) -> Molecule:
    """
    Generate a random molecule of type Molecule.
    """

    mol = Molecule()
    mol.atlist = generate_atom_list(
        config_generate,
        verbosity,
    )
    mol.num_atoms = np.sum(mol.atlist)
    mol.xyz, mol.ati = generate_coordinates(
        at=mol.atlist,
        scaling=config_generate.init_coord_scaling,
        dist_threshold=config_generate.dist_threshold,
        inc_scaling_factor=config_generate.increase_scaling_factor,
        verbosity=verbosity,
    )
    mol.charge = set_random_charge(mol.ati, verbosity)
    mol.uhf = 0
    mol.set_name_from_formula()

    if verbosity > 1:
        print(mol)

    return mol


def generate_atom_list(cfg: GenerateConfig, verbosity: int = 1) -> np.ndarray:
    """
    Generate a random molecule with a random number of atoms.
    """

    MAX_ELEM = 86
    # Define a new set of all elements that can be included
    set_all_elem = set(range(0, MAX_ELEM))
    if cfg.forbidden_elements:
        valid_elems = set_all_elem - set(cfg.forbidden_elements)
    else:
        valid_elems = set_all_elem

    natoms = np.zeros(
        102, dtype=int
    )  # 102 is the number of accessible elements in the periodic table

    # Reasoning for the parameters in the following sections:
    # - The number of the atoms added by default (DefaultRandom + AddOrganicAtoms)
    #   should match the minimum number of atoms in the molecule (if defined).
    # - The number of the atoms added by default (DefaultRandom + AddOrganicAtoms)
    #   should not exceed the maximum number of atoms in the molecule (if defined).
    # if the maximum number of atoms is not defined, the number of atoms added by default
    # should not exceed the minimum number of atoms + 20.
    # The current default value are the overall minimum if no other values are defined.
    # In general, the ratio of "default_random" atoms vs. "add_organic_atoms" atoms is 2:1.

    # With both 'add_random' and 'add_organic', we want to have ca. 60 % of the min_num_atoms
    # Reasoning: For an example of 6 atoms after 'add_random' and 'add_organic',
    #            the mean number of atoms added by 'add_hydrogen' is ca. 4
    # Thus: Round to the nearest integer of 0.6 * 2/3 = 0.4 of cfg.min_num_atoms
    if cfg.min_num_atoms:
        low_lim_default_random = round(cfg.min_num_atoms * (0.4))
        lim_organic = round(low_lim_default_random * 0.5)
    else:
        low_lim_default_random = 1  # default value if nothing is defined
        lim_organic = 5
    if cfg.max_num_atoms:
        max_lim_default_random = round((cfg.max_num_atoms + 10) * (0.4))
    else:
        if cfg.min_num_atoms:
            max_lim_default_random = round((cfg.min_num_atoms + 20) * (0.4))
        else:
            max_lim_default_random = 7  # default value if nothing is defined

    def add_random(min_adds: int, max_adds: int, min_nat: int, max_nat: int):
        """
        Default random atom generation.
        """
        numatoms_all = np.random.randint(
            min_adds, max_adds
        )  # with range(1, 7) -> mean value: 3.5
        for _ in range(numatoms_all):
            # Define the atom type to be added via a random choice from the set of valid elements
            ati = np.random.choice(list(valid_elems))
            if verbosity > 1:
                print(f"Adding atom type {ati}...")
            # Add a random number of atoms of the defined type
            natoms[ati] = natoms[ati] + np.random.randint(
                min_nat, max_nat
            )  # with range(0, 3) -> mean value: 1
            # max value of this section with commented settings: 12

    def add_organic(num_adds: int, min_nat: int, max_nat: int):
        """
        Add organic elements.
        """
        # Add Elements between B and F (5-9)
        for _ in range(num_adds):  # with range(5) -> mean value 1.5
            ati = np.random.randint(4, 10)
            if verbosity > 1:
                print(f"Adding atom type {ati}...")
            natoms[ati] = natoms[ati] + np.random.randint(
                min_nat, max_nat
            )  # with range(0, 3) -> mean value: 1
            # max value of this section with commented settings: 8

    def remove_group_onetwo():
        # > If too many alkaline and alkine earth metals are included, restart generation
        group_one_two = get_alkali_metals() + get_alkaline_earth_metals()
        nmetals = 0
        for i in group_one_two:
            nmetals += natoms[i]
        # reduce number of metals starting from 2, going to 55
        while nmetals > 3:
            for i in group_one_two:
                if natoms[i] > 0:
                    natoms[i] = natoms[i] - 1
                    if verbosity > 1:
                        print(f"Removing group 1/2 metal of type: {i}...")
                    nmetals -= 1
                if nmetals <= 3:
                    break

    def remove_metals():
        # If the sum of all other metals is larger than three, reduce the number of metals
        other_metals = (
            get_three_d_metals()
            + get_four_d_metals()
            + get_five_d_metals()
            + get_lanthanides()
        )
        n_othermetals = 0
        for i in other_metals:
            n_othermetals += natoms[i]
        while n_othermetals > 3:
            for i in other_metals:
                if natoms[i] > 0:
                    natoms[i] = natoms[i] - 1
                    if verbosity > 1:
                        print(f"Removing transition metal/LN of type: {i}...")
                    n_othermetals -= 1
                if n_othermetals <= 3:
                    break

    def add_hydrogen():
        # If no H is included, add H atoms
        if natoms[0] == 0:
            nat = np.sum(natoms)
            randint = np.random.rand()
            j = 1 + round(randint * nat * 1.2)
            natoms[0] = natoms[0] + j
            # Example: For 5 atoms at this point,
            # the mean number of added H atoms is (mean(1, 2, 3, 4, 5, 6))=3.5

    def check_min_max_atoms():
        # If the number of atoms is smaller than the minimum number of atoms, add atoms
        while np.sum(natoms) < cfg.min_num_atoms:
            if verbosity > 1:
                print(
                    f"Minimal number of atoms: {cfg.min_num_atoms}; Actual number of atoms: {np.sum(natoms)}.\nAdding atoms..."
                )
            ati = np.random.choice(list(valid_elems))
            max_limit = cfg.element_composition.get(ati, (None, None))[1]
            if max_limit is not None and natoms[ati] >= max_limit:
                continue
            natoms[ati] = natoms[ati] + 1
        # If the number of atoms is larger than the maximum number of atoms, remove atoms randomly
        tmp_count = 0
        while np.sum(natoms) > cfg.max_num_atoms:
            tmp_count += 1
            if tmp_count > 100:
                raise RuntimeError(
                    "Could not generate a molecule with the given constraints."
                )
            if verbosity > 1:
                print(
                    f"Max number of atoms: {cfg.max_num_atoms}; Actual number of atoms: {np.sum(natoms)}.\nRemoving atoms..."
                )
            # generate a list of all atom types that are included in the molecule with at least one atom
            # if the occurrence is > 1, add it multiple times to the list
            atom_list = []
            for i, count in enumerate(natoms):
                if count > 0:
                    atom_list.extend([i] * count)
            # randomly select an atom type from the list, thereby weighting the selection for reduction by the current occurrence
            # generate a random number between 0 and the number of atoms in the list
            random_index = np.random.randint(len(atom_list))
            i = atom_list[int(random_index)]
            if natoms[i] > 0:
                min_limit = cfg.element_composition.get(i, (None, None))[0]
                if verbosity > 1:
                    print(f"Trying to remove atom type {i}...")
                if min_limit is None or natoms[i] > min_limit:
                    natoms[i] = natoms[i] - 1

    def check_composition():
        # Align with the given element_composition:
        # CAUTION: The setting to min/max count may violate the metal count restrictions
        for elem, count_range in cfg.element_composition.items():
            min_count, max_count = count_range
            if min_count is not None and natoms[elem] < min_count:
                natoms[elem] = min_count
            elif max_count is not None and natoms[elem] > max_count:
                natoms[elem] = max_count

    ### ACTUAL WORKFLOW START ###
    # Add a random number of atoms of random types
    add_random(low_lim_default_random, max_lim_default_random, 0, 3)
    # Check for too many group 1 and 2 metals
    remove_group_onetwo()
    # Check for too many transition and lanthanide metals
    remove_metals()
    # Add organic elements (B, C, N, O, F)
    add_organic(lim_organic, 0, 3)
    # Add hydrogen if not included
    add_hydrogen()
    # Check if pre-defined atom type counts are within the defined limits
    check_composition()
    # Check if the number of atoms is within the defined limits
    check_min_max_atoms()
    ### ACTUAL WORKFLOW END ###

    return natoms


def generate_coordinates(
    at: np.ndarray,
    scaling: float,
    dist_threshold: float,
    inc_scaling_factor: float = 1.3,
    verbosity: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate random coordinates for a molecule.
    """

    # eff_scaling is a deep copy of scaling
    eff_scaling = copy.deepcopy(scaling)
    xyz, ati = generate_random_coordinates(at)
    xyz = xyz * eff_scaling
    # do while check_distances is False
    while not check_distances(xyz, dist_threshold):
        if verbosity > 1:
            print(
                f"Distance check failed. Increasing expansion factor by {inc_scaling_factor}..."
            )
        xyz, ati = generate_random_coordinates(at)
        eff_scaling = eff_scaling * inc_scaling_factor
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


def get_alkali_metals() -> list[int]:
    """
    Get the atomic numbers of alkali metals.
    """
    alkali = [2, 10, 18, 36, 54]
    return alkali


def get_alkaline_earth_metals() -> list[int]:
    """
    Get the atomic numbers of alkaline earth metals.
    """
    alkaline = [3, 11, 19, 37, 55]
    return alkaline


def get_three_d_metals() -> list[int]:
    """
    Get the atomic numbers of three d metals.
    """
    threedmetals = list(range(20, 30))

    return threedmetals


def get_four_d_metals() -> list[int]:
    """
    Get the atomic numbers of four d metals.
    """

    fourdmetals = list(range(38, 48))
    return fourdmetals


def get_five_d_metals() -> list[int]:
    """
    Get the atomic numbers of five d metals.
    """
    fivedmetals = list(range(71, 80))
    return fivedmetals


def get_lanthanides() -> list[int]:
    """
    Get the atomic numbers of lanthanides.
    """
    lanthanides = list(range(56, 71))
    return lanthanides
