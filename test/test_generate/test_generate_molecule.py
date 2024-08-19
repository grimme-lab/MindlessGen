import pytest
import numpy as np
from mindlessgen.molecules import (  # type: ignore
    generate_atom_list,
    get_alkali_metals,
    get_alkaline_earth_metals,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
    Molecule,
)
from mindlessgen.prog import GenerateConfig  # type: ignore


@pytest.fixture
def default_generate_config():
    """Fixture to provide a default GenerateConfig object."""
    config = GenerateConfig()
    config.min_num_atoms = 5
    config.max_num_atoms = 15
    config.init_coord_scaling = 1.0
    config.dist_threshold = 0.8
    config.increase_scaling_factor = 1.2
    return config


@pytest.fixture
def default_molecule():
    """Fixture to provide a default Molecule object."""
    return Molecule()


@pytest.mark.parametrize("min_atoms,max_atoms", [(5, 15), (10, 20), (2, 10)])
def test_generate_atom_list_combined(min_atoms, max_atoms, default_generate_config):
    """Test the generate_atom_list function with various atom ranges and validate metal constraints."""
    # Update the generate config with the parametrized values
    default_generate_config.min_num_atoms = min_atoms
    default_generate_config.max_num_atoms = max_atoms

    # Generate the atom list
    atom_list = generate_atom_list(default_generate_config, verbosity=1)

    # Common assertions for atom count range
    assert isinstance(atom_list, np.ndarray)
    assert np.sum(atom_list) >= min_atoms
    assert np.sum(atom_list) <= max_atoms

    # Additional assertions from the second test
    assert atom_list.shape == (102,)
    assert np.sum(atom_list) > 0

    # Check that for the transition and lanthanide metals, the occurrence is never greater than 3
    all_metals = (
        get_three_d_metals()
        + get_four_d_metals()
        + get_five_d_metals()
        + get_lanthanides()
    )
    for z in all_metals:
        assert atom_list[z] <= 3

    # Check that the sum of alkali and alkaline earth metals is never greater than 3
    alkmetals = get_alkali_metals() + get_alkaline_earth_metals()
    assert np.sum([atom_list[z] for z in alkmetals]) <= 3


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
