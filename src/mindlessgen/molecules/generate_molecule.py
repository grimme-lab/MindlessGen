"""
This module generates a random molecule with a random number of atoms.
"""

import copy
import numpy as np
from ..prog import GenerateConfig
from .molecule import Molecule
from .refinement import get_cov_radii, COV_RADII
from .miscellaneous import (
    set_random_charge,
    calculate_protons,
    get_alkali_metals,
    get_alkaline_earth_metals,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
    get_actinides,
    calculate_ligand_electrons,
    calculate_uhf,
)


MAX_ELEM = 86


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
        inc_scaling_factor=config_generate.increase_scaling_factor,
        verbosity=verbosity,
        scale_minimal_distance=config_generate.scale_minimal_distance,
    )
    if config_generate.contract_coords:
        mol.xyz = contract_coordinates(
            xyz=mol.xyz,
            ati=mol.ati,
            scale_minimal_distance=config_generate.scale_minimal_distance,
        )
    if config_generate.molecular_charge is not None:
        mol.charge = config_generate.molecular_charge
        mol.uhf = calculate_uhf(mol.atlist)
    else:
        mol.charge, mol.uhf = set_random_charge(mol.ati, verbosity)
    mol.set_name_from_formula()
    if verbosity > 1:
        print(mol)

    return mol


# Taken from mlmgen Fortran project.
# TODO: Make this procedure more object-oriented:
#       Create a generator base class that contains some basic functions.
#       Employ more specific generator classes that have superseeding functions
#       whenever needed.
def generate_atom_list(cfg: GenerateConfig, verbosity: int = 1) -> np.ndarray:
    """
    Generate a random molecule with a random number of atoms.
    """
    # initialize a default random number generator
    rng = np.random.default_rng()

    # Define a new set of all elements that can be included
    set_all_elem = set(range(0, MAX_ELEM))
    if cfg.forbidden_elements:
        valid_elems = set_all_elem - set(cfg.forbidden_elements)
    else:
        valid_elems = set_all_elem

    natoms = np.zeros(103, dtype=int)  # Support for up to element 103 (Lr)

    if cfg.fixed_composition:
        # Set the number of atoms for the fixed composition
        for elem, count_range in cfg.element_composition.items():
            natoms[elem] = count_range[0]
        if verbosity > 1:
            print(
                "Setting the number of atoms for the fixed composition. "
                + f"Returning: \n{natoms}\n"
            )
        # If the molecular charge is defined, and a fixed element composition is defined, check if the electrons are even. If not raise an error.
        if cfg.molecular_charge:
            protons = calculate_protons(natoms)
            nel = protons - cfg.molecular_charge
            f_elem = any(
                count > 0 and (i in get_lanthanides() or i in get_actinides())
                for i, count in enumerate(natoms)
            )
            if (f_elem and calculate_ligand_electrons(natoms, nel) % 2 != 0) or (
                not f_elem and nel % 2 != 0
            ):
                raise ValueError(
                    "Both fixed charge and fixed composition are defined. "
                    + "Please only define one of them."
                    + "Or ensure that the fixed composition is closed shell."
                )
        return natoms

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
        numatoms_all = rng.integers(
            low=min_adds, high=max_adds
        )  # with range(1, 7) -> mean value: 3.5
        for _ in range(numatoms_all):
            # Define the atom type to be added via a random choice from the set of valid elements
            ati = rng.choice(list(valid_elems))
            if verbosity > 1:
                print(f"Adding atom type {ati}...")
            # Add a random number of atoms of the defined type
            natoms[ati] = natoms[ati] + rng.integers(
                low=min_nat, high=max_nat
            )  # with range(0, 3) -> mean value: 1
            # max value of this section with commented settings: 12

    def add_organic(num_adds: int, min_nat: int, max_nat: int) -> None:
        """
        Add organic elements.
        """
        # Add Elements between B and F (5-9)
        valid_organic: list[int] = []
        for organic_index in range(4, 10):
            if organic_index in valid_elems:
                valid_organic.append(organic_index)
        if not valid_organic:
            return
        for _ in range(num_adds):  # with range(5) -> mean value 1.5
            # go through the elements B to F (4-9 in 0-based indexing)
            ati = rng.choice(valid_organic)
            if verbosity > 1:
                print(f"Adding atom type {ati}...")
            natoms[ati] = natoms[ati] + rng.integers(
                low=min_nat, high=max_nat
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
            randint = rng.random()
            j = 1 + round(randint * nat * 1.2)
            natoms[0] = natoms[0] + j
            # Example: For 5 atoms at this point,
            # the mean number of added H atoms is (mean(1, 2, 3, 4, 5, 6))=3.5

    def check_min_max_atoms():
        # If the number of atoms is smaller than the minimum number of atoms, add atoms
        while np.sum(natoms) < cfg.min_num_atoms:
            if verbosity > 1:
                print(
                    f"Minimal number of atoms: {cfg.min_num_atoms}; "
                    + f"Actual number of atoms: {np.sum(natoms)}.\nAdding atoms..."
                )
            ati = rng.choice(list(valid_elems))
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
                    f"Max number of atoms: {cfg.max_num_atoms}; "
                    + f"Actual number of atoms: {np.sum(natoms)}.\nRemoving atoms..."
                )
            # generate a list of all atom types that are included in the molecule
            # with at least one atom
            # if the occurrence is > 1, add it multiple times to the list
            atom_list = []
            for i, count in enumerate(natoms):
                if count > 0:
                    atom_list.extend([i] * count)
            # randomly select an atom type from the list, thereby weighting the selection
            # for reduction by the current occurrence
            # generate a random number between 0 and the number of atoms in the list
            random_index = rng.integers(0, len(atom_list))
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
            # define the random number of atoms to be added
            if min_count is None:
                min_count = 0
            if max_count is None:
                max_count = cfg.max_num_atoms
                # 50 % of the maximally allowed number of atoms
            added_atoms_from_composition = rng.integers(
                low=min_count, high=max_count, endpoint=True
            )
            natoms[elem] = added_atoms_from_composition

    ### ACTUAL WORKFLOW START ###
    # Add a random number of atoms of random types
    add_random(low_lim_default_random, max_lim_default_random, 0, 3)
    # Check for too many group 1 and 2 metals
    remove_group_onetwo()
    # Check for too many transition and lanthanide metals
    remove_metals()
    # Add organic elements (B, C, N, O, F)
    add_organic(num_adds=lim_organic, min_nat=0, max_nat=3)
    # Add hydrogen if not included
    # execute only if hydrogen is included in the valid elements
    if 0 in valid_elems:
        add_hydrogen()
    # Check if pre-defined atom type counts are within the defined limits
    check_composition()
    # Check if the number of atoms is within the defined limits
    check_min_max_atoms()
    # If the molecule is not closed shell, add an atom to ensure a closed shell system
    if cfg.molecular_charge is not None:
        protons = calculate_protons(natoms)
        nel = protons - cfg.molecular_charge
        natoms = fixed_charge_correction(cfg, natoms, nel, valid_elems, verbosity)
    ### ACTUAL WORKFLOW END ###

    return natoms


def generate_coordinates(
    at: np.ndarray,
    scaling: float,
    inc_scaling_factor: float = 1.3,
    verbosity: int = 1,
    scale_minimal_distance: float = 0.8,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate random coordinates for a molecule.
    """

    # eff_scaling is a deep copy of scaling
    eff_scaling = copy.deepcopy(scaling)
    xyz, ati = generate_random_coordinates(at)
    xyz = xyz * eff_scaling
    # do while check_distances is False
    while not check_distances(xyz, ati, scale_minimal_distance=scale_minimal_distance):
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
    rng = np.random.default_rng()
    for elem, count in enumerate(at):
        for m in range(count):
            # different rules for hydrogen
            if elem == 0:
                xyz[numatoms + m, :] = rng.random(3) * 3 - 1.5
            else:
                xyz[numatoms + m, :] = rng.random(3) * 2 - 1

            atilist.append(elem)

        numatoms += count

    ati = np.array(atilist, dtype=int)

    return xyz, ati


def contract_coordinates(
    xyz: np.ndarray, ati: np.ndarray, scale_minimal_distance: float
) -> np.ndarray:
    """
    Pull the atoms towards the origin.
    """
    # Initialize the old coordinates as an array of zeros
    xyz_old: np.ndarray = np.zeros_like(xyz)
    cycle = 0
    # Break if the coordinates do not change
    while not np.array_equal(xyz_old, xyz):
        cycle += 1
        if cycle > 2500:
            raise RuntimeError(
                "Could not contract the coordinates in a reasonable amount of cycles."
            )
        xyz_old = xyz.copy()
        # Go through the atoms dimension of the xyz array in a reversed order
        # Justification: First atoms are most likely hydrogen atoms, which should be moved last
        for i in range(len(xyz) - 1, -1, -1):
            atom_xyz = xyz[i]
            atom_xyz_norm = np.linalg.norm(atom_xyz)
            normalized_atom_xyz = atom_xyz / atom_xyz_norm
            # Shift the atom only if it is closer to the origin after the shift
            if atom_xyz_norm > 0.1:
                # Pull the atom towards the origin
                xyz[i] -= normalized_atom_xyz * 0.2
                # When the check_distances function returns False, reset the atom coordinates
                if not check_distances(xyz, ati, scale_minimal_distance):
                    xyz[i] = xyz_old[i]
    return xyz


def check_distances(
    xyz: np.ndarray, ati: np.ndarray, scale_minimal_distance: float
) -> bool:
    """
    Check if the distances between atoms are larger than a threshold.
    """
    # go through the atoms dimension of the xyz array
    for i in range(xyz.shape[0] - 1):
        for j in range(i + 1, xyz.shape[0]):
            r = np.linalg.norm(xyz[i, :] - xyz[j, :])
            sum_radii = get_cov_radii(ati[i], COV_RADII) + get_cov_radii(
                ati[j], COV_RADII
            )
            if r < scale_minimal_distance * sum_radii:
                return False
    return True


def fixed_charge_correction(
    cfg: GenerateConfig,
    natoms: np.ndarray,
    nel: int,
    valid_elems: set[int],
    verbosity: int,
) -> np.ndarray:
    """
    Correct the number of electrons if a fixed charge is given and the molecule is not closed shell.
    """
    f_elem = any(
        count > 0 and (i in get_lanthanides() or i in get_actinides())
        for i, count in enumerate(natoms)
    )
    if f_elem:
        ligand_electrons = calculate_ligand_electrons(natoms, nel)
        # If f block elements are included, correct only if the remaning ligand protons are uneven
        if ligand_electrons % 2 != 0:
            natoms = fixed_charge_elem_correction(cfg, natoms, valid_elems, verbosity)
            return natoms
    # If f block elements are not included, correct if the number of electrons is uneven
    elif nel % 2 != 0:
        natoms = fixed_charge_elem_correction(cfg, natoms, valid_elems, verbosity)
        return natoms
    return natoms


def fixed_charge_elem_correction(
    cfg: GenerateConfig,
    natoms: np.ndarray,
    valid_elems: set[int],
    verbosity: int,
) -> np.ndarray:
    """
    Correct the number of atoms if the number of electrons is odd and a molecular charge is set.
    """
    num_atoms = np.sum(natoms)
    # All other elements
    rng = np.random.default_rng()
    odd_atoms = np.array([elem for elem in valid_elems if elem % 2 == 0], dtype=int)
    random_odd_atoms = rng.permutation(odd_atoms)
    if 0 in valid_elems:
        random_odd_atoms = np.insert(random_odd_atoms, 0, 0)
    for random_elem in random_odd_atoms:
        min_count, max_count = cfg.element_composition.get(
            random_elem, (0, cfg.max_num_atoms)
        )
        if min_count is None:
            min_count = 0
        if max_count is None:
            max_count = cfg.max_num_atoms
        # Check if adding or removing the random element is possible
        if natoms[random_elem] < max_count and num_atoms < cfg.max_num_atoms:
            natoms[random_elem] += 1
            if verbosity > 1:
                print(f"Adding atom type {random_elem} for charge...")
            return natoms
        if natoms[random_elem] > min_count and num_atoms > cfg.min_num_atoms:
            natoms[random_elem] -= 1
            if verbosity > 1:
                print(f"Removing atom type {random_elem} for charge...")
            return natoms
    raise RuntimeError("Could not correct the odd number of electrons.")
