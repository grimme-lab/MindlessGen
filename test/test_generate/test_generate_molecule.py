import pytest
import numpy as np
from mindlessgen.molecules import (  # type: ignore
    generate_atom_list,
    generate_random_molecule,
    generate_coordinates,
    check_distances,
    get_alkali_metals,
    get_alkaline_earth_metals,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
)
from mindlessgen.prog import GenerateConfig, ConfigManager  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore


@pytest.fixture
def default_generate_config():
    """Fixture to provide a default GenerateConfig object."""
    config = GenerateConfig()
    return config


@pytest.mark.parametrize("min_atoms,max_atoms", [(5, 15), (10, 20), (2, 10)])
def test_generate_atom_list(min_atoms, max_atoms, default_generate_config):
    """Test the generate_atom_list function with various atom ranges and validate metal constraints."""
    # Update the generate config with the parametrized values
    default_generate_config.min_num_atoms = min_atoms
    default_generate_config.max_num_atoms = max_atoms

    # Generate the atom list
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Common assertions for atom count range
    assert isinstance(atom_list, np.ndarray)
    assert atom_list.shape == (102,)
    assert np.sum(atom_list) > 0
    assert np.sum(atom_list) >= min_atoms
    assert np.sum(atom_list) <= max_atoms

    # Check that the sum of transition and lanthanide metals is never greater than 3
    all_metals = (
        get_three_d_metals()
        + get_four_d_metals()
        + get_five_d_metals()
        + get_lanthanides()
    )
    assert np.sum([atom_list[z] for z in all_metals]) <= 3

    # Check that the sum of alkali and alkaline earth metals is never greater than 3
    alkmetals = get_alkali_metals() + get_alkaline_earth_metals()
    assert np.sum([atom_list[z] for z in alkmetals]) <= 3


# Test the element composition property of the GenerateConfig class
def test_generate_config_element_composition(default_generate_config):
    """Test the element composition property of the GenerateConfig class."""
    default_generate_config.min_num_atoms = 10
    default_generate_config.max_num_atoms = 15
    default_generate_config.element_composition = "C:2-2, N:3-3, O:1-1"
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Check that the atom list contains the correct number of atoms for each element
    assert atom_list[5] == 2
    assert atom_list[6] == 3
    assert atom_list[7] == 1


# Test the forbidden_elements property of the GenerateConfig class
def test_generate_config_forbidden_elements(default_generate_config):
    """Test the forbidden_elements property of the GenerateConfig class."""
    default_generate_config.min_num_atoms = 10
    default_generate_config.max_num_atoms = 15
    default_generate_config.forbidden_elements = "7-9, 19-*"
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Check that the atom list does not contain any forbidden elements
    assert np.sum([atom_list[z] for z in range(6, 9)]) == 0
    assert np.sum([atom_list[z] for z in range(18, 102)]) == 0


@pytest.mark.parametrize(
    "min_num_atoms, max_num_atoms, element_composition, expected_error, expected_atom_counts",
    [
        # Case 1: Summed min atoms larger than max atoms
        (10, 10, "C:5-10, H:6-10", ValueError, None),
        # Case 2: Fixed composition outside of defined limits
        (5, 10, "C:3-3, H:3-3, O:5-5", ValueError, None),
        # Case 3: Fixed composition within defined limits
        (
            5,
            10,
            "C:3-3, H:2-2",
            None,
            {5: 3, 0: 2},  # Carbon (C) and Hydrogen (H)
        ),
        # Case 4: Not all elements in composition are fixed
        (
            5,
            10,
            "C:2-4, H:2-2",
            None,
            {5: (2, 4), 0: 2},  # Carbon (C) in range, Hydrogen (H) fixed
        ),
    ],
)
def test_generate_atom_list_with_composition(
    min_num_atoms,
    max_num_atoms,
    element_composition,
    expected_error,
    expected_atom_counts,
):
    """
    Parametrized test for checking element composition conditions.
    """
    config = GenerateConfig()
    config.min_num_atoms = min_num_atoms
    config.max_num_atoms = max_num_atoms
    config.element_composition = element_composition

    if expected_error:
        with pytest.raises(expected_error):
            generate_atom_list(config, verbosity=1)
    else:
        atom_list = generate_atom_list(config, verbosity=1)
        for elem, count in expected_atom_counts.items():
            if isinstance(count, tuple):
                assert count[0] <= atom_list[elem] <= count[1]
            else:
                assert atom_list[elem] == count


def test_get_alkali_metals():
    """Test the get_alkali_metals function."""
    alkali_metals = get_alkali_metals()
    assert isinstance(alkali_metals, list)
    assert all(isinstance(e, int) for e in alkali_metals)


def test_get_alkaline_earth_metals():
    """Test the get_alkaline_earth_metals function."""
    alkaline_earth_metals = get_alkaline_earth_metals()
    assert isinstance(alkaline_earth_metals, list)
    assert all(isinstance(e, int) for e in alkaline_earth_metals)


def test_get_three_d_metals():
    """Test the get_three_d_metals function."""
    three_d_metals = get_three_d_metals()
    assert isinstance(three_d_metals, list)
    assert all(isinstance(e, int) for e in three_d_metals)


def test_get_four_d_metals():
    """Test the get_four_d_metals function."""
    four_d_metals = get_four_d_metals()
    assert isinstance(four_d_metals, list)
    assert all(isinstance(e, int) for e in four_d_metals)


def test_get_five_d_metals():
    """Test the get_five_d_metals function."""
    five_d_metals = get_five_d_metals()
    assert isinstance(five_d_metals, list)
    assert all(isinstance(e, int) for e in five_d_metals)


def test_get_lanthanides():
    """Test the get_lanthanides function."""
    lanthanides = get_lanthanides()
    assert isinstance(lanthanides, list)
    assert all(isinstance(e, int) for e in lanthanides)


def test_generate_molecule() -> None:
    """
    Test the generation of an array of atomic numbers.
    """
    # create a ConfigManager object with verbosity set to 0
    config = ConfigManager()
    config.generate.forbidden_elements = "57-71"
    config.general.verbosity = 0
    mol = generate_random_molecule(config.generate, config.general.verbosity)

    assert mol.num_atoms > 0
    assert mol.num_atoms == np.sum(mol.atlist)
    assert mol.num_atoms == len(mol.xyz)
    assert mol.num_atoms == len(mol.ati)
    assert mol.uhf == 0
    assert mol.charge is not None
    # assert that sum of absolute value of mol.xyz is greater than 0
    assert np.sum(np.abs(mol.xyz)) > 0


def test_generate_coordinates() -> None:
    """
    Test the generation of coordinates.
    """
    mol = Molecule()
    # create an empty array with dimension 1 and length 86
    mol.atlist = np.zeros(102, dtype=int)
    assert mol.atlist.shape == (102,)

    # set the first element to 4 and the sixth element to 2
    mol.atlist[0] = 4
    mol.atlist[5] = 2

    mol.xyz, mol.ati = generate_coordinates(mol.atlist, 3.0, 1.2)
    # assert that the shape of mol.xyz is (6, 3)
    assert mol.xyz.shape == (6, 3)
    assert mol.num_atoms == 6
    assert mol.num_atoms == np.sum(mol.atlist)
    assert mol.num_atoms == len(mol.xyz)
    assert mol.num_atoms == len(mol.ati)


@pytest.mark.parametrize(
    "xyz, threshold, expected, description",
    [
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            0.5,
            True,
            "Two atoms with distance greater than threshold (1.0 > 0.5)",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]]),
            0.5,
            False,
            "Two atoms with distance less than threshold (0.4 < 0.5)",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            0.5,
            True,
            "Three atoms in a line with distances greater than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            0.5,
            False,
            "Three atoms with one pair close together: distance between first two is less than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0]]),
            0.5,
            True,
            "Single atom, no distances to compare",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            0.5,
            False,
            "Two atoms at identical positions: distance is zero, less than threshold",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            1.7320,
            True,
            "Two atoms with diagonal distance just above threshold (sqrt(3) ≈ 1.732)",
        ),
    ],
    ids=[
        "far_apart",
        "close_together",
        "three_in_line",
        "three_with_one_close",
        "single_atom",
        "two_identical",
        "diagonal_distance",
    ],
)
def test_check_distances(xyz, threshold, expected, description):
    assert check_distances(xyz, threshold) == expected


def test_generate_atom_list_min_larger_than_max(default_generate_config):
    """Test the generate_atom_list function when the minimum number of atoms is larger than the maximum number of atoms."""
    default_generate_config.min_num_atoms = 10
    default_generate_config.max_num_atoms = 5

    with pytest.raises(ValueError):
        generate_atom_list(default_generate_config, verbosity=1)


# Test to ensure non-integer values for min/max atoms raise errors
@pytest.mark.parametrize(
    "min_atoms, max_atoms, expected_error",
    [
        ("five", 10, TypeError),  # Non-integer min atoms
        (5, "ten", TypeError),  # Non-integer max atoms
        (2.5, 15, TypeError),  # Float as min atoms
        (5, 20.5, TypeError),  # Float as max atoms
    ],
)
def test_invalid_min_max_atoms(
    min_atoms, max_atoms, expected_error, default_generate_config
):
    """Test for non-integer min/max atom values raising TypeErrors."""
    with pytest.raises(expected_error):
        default_generate_config.min_num_atoms = min_atoms
        default_generate_config.max_num_atoms = max_atoms


# Edge case where forbidden elements overlap allowed range
@pytest.mark.parametrize(
    "forbidden_elements, expected_atoms",
    [
        ("1-5", [0, 1, 2, 3, 4]),  # Banned elements within the default organic range
        ("1, 6, 7", [0, 5, 6]),  # Specific elements banned
    ],
)
def test_generate_atom_list_with_overlapping_forbidden_elements(
    forbidden_elements, expected_atoms, default_generate_config
):
    """Test generate_atom_list when forbidden elements overlap with allowed ranges."""
    default_generate_config.forbidden_elements = forbidden_elements
    default_generate_config.min_num_atoms = 5
    default_generate_config.max_num_atoms = 15
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Ensure forbidden elements are not present in the atom list
    assert np.sum([atom_list[z] for z in expected_atoms]) == 0


# Test behavior when composition is empty but min/max are set
def test_generate_atom_list_with_empty_composition(default_generate_config):
    """Ensure empty compositions don't lead to unexpected behaviors."""
    default_generate_config.element_composition = ""
    default_generate_config.min_num_atoms = 5
    default_generate_config.max_num_atoms = 10
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Ensure some atoms are still generated within min and max limits
    assert np.sum(atom_list) >= default_generate_config.min_num_atoms
    assert np.sum(atom_list) <= default_generate_config.max_num_atoms


# Test for element compositions with zero ranges
def test_generate_atom_list_zero_composition(default_generate_config):
    """Test generate_atom_list when compositions have zero counts."""
    default_generate_config.element_composition = "C:0-0, N:0-0, O:0-0"
    default_generate_config.min_num_atoms = 5
    default_generate_config.max_num_atoms = 10
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Ensure atoms in these ranges are indeed zero
    assert atom_list[5] == 0  # C
    assert atom_list[6] == 0  # N
    assert atom_list[7] == 0  # O


# Check hydrogen addition when it should/shouldn't occur
@pytest.mark.parametrize(
    "forbidden_elements, should_contain_hydrogen",
    [
        ("1", False),  # Hydrogen forbidden
        ("", True),  # No forbidden elements
    ],
)
def test_hydrogen_addition(
    forbidden_elements, should_contain_hydrogen, default_generate_config
):
    """Test hydrogen addition based on forbidden elements."""
    default_generate_config.forbidden_elements = forbidden_elements
    default_generate_config.min_num_atoms = 5
    default_generate_config.max_num_atoms = 15
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    if should_contain_hydrogen:
        assert atom_list[0] > 0
    else:
        np.testing.assert_equal(atom_list[0], 0)
