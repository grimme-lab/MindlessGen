import pytest
from mindlessgen.prog import (  # type: ignore
    GeneralConfig,
    GenerateConfig,
    RefineConfig,
    XTBConfig,
    ORCAConfig,
    PostProcessConfig,
    TURBOMOLEConfig,
    SymmetrizationConfig,
)


# Tests for GeneralConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("verbosity", 2, "high", TypeError),
        ("verbosity", 2, 4, ValueError),
        ("max_cycles", 50, -1, ValueError),
        ("max_cycles", 50, "100", TypeError),
        ("print_config", True, "yes", TypeError),
        ("parallel", 4, 0, ValueError),
        ("parallel", 4, "four", TypeError),
        ("num_molecules", 2, 0, ValueError),
        ("num_molecules", 2, "two", TypeError),
        ("postprocess", True, "yes", TypeError),
        ("write_xyz", True, "yes", TypeError),
        ("symmetrization", False, "no", TypeError),
    ],
)
def test_general_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = GeneralConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("verbosity", 1),
        ("max_cycles", 200),
        ("print_config", False),
        ("parallel", 4),
        ("num_molecules", 1),
        ("postprocess", False),
        ("write_xyz", True),
        ("symmetrization", False),
    ],
)
def test_general_config_default_values(property_name, initial_value):
    """
    Test default values for GeneralConfig properties.
    """
    config = GeneralConfig()
    assert getattr(config, property_name) == initial_value


# Tests for GenerateConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("min_num_atoms", 5, 0, ValueError),
        ("min_num_atoms", 5, "two", TypeError),
        ("max_num_atoms", 80, -10, ValueError),
        ("max_num_atoms", 80, None, TypeError),
        ("init_coord_scaling", 1.0, -0.5, ValueError),
        ("init_coord_scaling", 1.0, "1.0", TypeError),
        ("increase_scaling_factor", 1.1, 0.0, ValueError),
        ("increase_scaling_factor", 1.1, "1.1", TypeError),
        ("scale_fragment_detection", 1.25, "1.25", TypeError),
        ("scale_fragment_detection", 1.25, 0.0, ValueError),
        ("scale_minimal_distance", 0.8, "0.8", TypeError),
        ("scale_minimal_distance", 0.8, -1.0, ValueError),
        ("contract_coords", True, "true", TypeError),
        ("molecular_charge", 1, [], TypeError),
        ("molecular_charge", 2, "not_a_number", ValueError),
        ("fixed_composition", False, "no", TypeError),
    ],
)
def test_generate_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    """
    Test property setters for GenerateConfig class.
    """
    config = GenerateConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


# For the property: "element_composition" and "forbidden_elements", introduce special tests
# Check that setting element_composition with "C:1-10" leads to "{5: (1, 10)}"
# Check that setting element_composition with "C:1-10, H:2-*" leads to "{5: (1, 10), 0: (2, None)}"
@pytest.mark.parametrize(
    "property_name, valid_value, expected_value",
    [
        ("element_composition", "C:1-10", {5: (1, 10)}),
        ("element_composition", "C:1-10, H:2-*", {5: (1, 10), 0: (2, None)}),
        # additional test for giving the element_composition directly as a dictionary
        ("element_composition", {5: (1, 10), 0: (2, None)}, {5: (1, 10), 0: (2, None)}),
        ("forbidden_elements", "6,1", [0, 5]),
        (
            "forbidden_elements",
            "86-*",
            [85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102],
        ),
        # additional test for giving the forbidden_elements directly as a list
        ("forbidden_elements", [0, 5], [0, 5]),
    ],
)
def test_generate_config_element_composition(
    property_name, valid_value, expected_value
):
    config = GenerateConfig()
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == expected_value


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("min_num_atoms", 5),
        ("max_num_atoms", 10),
        ("init_coord_scaling", 3.0),
        ("increase_scaling_factor", 1.1),
        ("element_composition", {}),
        ("forbidden_elements", None),
        ("scale_fragment_detection", 1.25),
        ("scale_minimal_distance", 0.8),
        ("contract_coords", True),
        ("molecular_charge", None),
        ("fixed_composition", False),
    ],
)
def test_generate_config_default_values(property_name, initial_value):
    """
    Test default values for GenerateConfig properties.
    """
    config = GenerateConfig()
    assert getattr(config, property_name) == initial_value


# Test for fixed charge
@pytest.mark.parametrize(
    "property_name, valid_value, expected_value",
    [
        ("molecular_charge", 0, 0),
        ("molecular_charge", "1", 1),
        ("molecular_charge", "none", None),
        ("molecular_charge", "", None),
    ],
)
def test_generate_config_molecular_charge(property_name, valid_value, expected_value):
    config = GenerateConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == expected_value


@pytest.mark.parametrize(
    "property_name, invalid_value, expected_exception",
    [
        ("molecular_charge", "two", ValueError),
        ("molecular_charge", "1.0", ValueError),
    ],
)
def test_generate_config_molecular_charge_invalid(
    property_name, invalid_value, expected_exception
):
    config = GenerateConfig()

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


# Tests for RefineConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("max_frag_cycles", 200, -1, ValueError),
        ("max_frag_cycles", 200, "100", TypeError),
        ("engine", "xtb", 123, TypeError),
        ("engine", "xtb", "g16", ValueError),
        ("hlgap", 0.5, "0.5", TypeError),
        ("hlgap", 0.5, -1.0, ValueError),
        ("debug", False, "false", TypeError),
        ("ncores", 2, "two", TypeError),
    ],
)
def test_refine_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = RefineConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("max_frag_cycles", 10),
        ("engine", "xtb"),
        ("hlgap", 0.5),
        ("debug", False),
        ("ncores", 2),
    ],
)
def test_refine_config_default_values(property_name, initial_value):
    """
    Test default values for RefineConfig properties.
    """
    config = RefineConfig()
    assert getattr(config, property_name) == initial_value


# Tests for XTBConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("xtb_path", "path/to/xtb", 123, TypeError),
        ("xtb_path", "path/to/xtb", None, TypeError),
        ("level", 2, "high", TypeError),
        ("level", 2, -1, ValueError),
        ("level", 1, 3, ValueError),
    ],
)
def test_xtb_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    """
    Test property setters for XTBConfig class.
    """
    config = XTBConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("xtb_path", "xtb"),
        ("level", 2),
    ],
)
def test_xtb_config_default_values(property_name, initial_value):
    """
    Test default values for XTBConfig properties.
    """
    config = XTBConfig()
    assert getattr(config, property_name) == initial_value


# Tests for ORCAConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("orca_path", "path/to/orca", 123, TypeError),
        ("orca_path", "path/to/orca", None, TypeError),
        ("functional", "PBE0", 123, TypeError),
        ("basis", "def2-TZVP", False, TypeError),
        ("gridsize", 2, "fine", TypeError),
        ("gridsize", 2, 4, ValueError),
        ("scf_cycles", 150, "many", TypeError),
        ("scf_cycles", 150, 0, ValueError),
    ],
)
def test_orca_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = ORCAConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("orca_path", "orca"),
        ("functional", "PBE"),
        ("basis", "def2-SVP"),
        ("gridsize", 1),
        ("scf_cycles", 100),
    ],
)
def test_orca_config_default_values(property_name, initial_value):
    """
    Test default values for ORCAConfig properties.
    """
    config = ORCAConfig()
    assert getattr(config, property_name) == initial_value


# Tests for PostProcessConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("engine", "orca", 123, TypeError),
        ("engine", "orca", "g16", ValueError),
        ("optimize", True, "yes", TypeError),
        ("opt_cycles", 10, -1, ValueError),
        ("opt_cycles", 10, 1.5, TypeError),
        ("opt_cycles", "none", "maybe", ValueError),
        ("debug", False, "false", TypeError),
        ("ncores", 4, "four", TypeError),
    ],
)
def test_postprocess_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = PostProcessConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    expected = (
        None
        if (property_name == "opt_cycles" and valid_value == "none")
        else valid_value
    )
    assert getattr(config, property_name) == expected

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("engine", "orca"),
        ("opt_cycles", None),
        ("optimize", True),
        ("debug", False),
        ("ncores", 4),
    ],
)
def test_postprocess_config_default_values(property_name, initial_value):
    """
    Test default values for PostProcessConfig properties.
    """
    config = PostProcessConfig()
    assert getattr(config, property_name) == initial_value


@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("ridft_path", "ridft", 123, TypeError),
        ("jobex_path", "jobex", None, TypeError),
        ("functional", "pbe", 123, TypeError),
        ("basis", "def2-SVP", False, TypeError),
        ("scf_cycles", 150, "many", TypeError),
        ("scf_cycles", 150, 0, ValueError),
    ],
)
def test_turbomole_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = TURBOMOLEConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("ridft_path", "ridft"),
        ("jobex_path", "jobex"),
        ("functional", "pbe"),
        ("basis", "def2-SVP"),
        ("scf_cycles", 100),
    ],
)
def test_turbomole_config_default_values(property_name, initial_value):
    """
    Test default values for TURBOMOLEConfig properties.
    """
    config = TURBOMOLEConfig()
    assert getattr(config, property_name) == initial_value


@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("distance", 3.0, -1.0, ValueError),
        ("distance", 3.0, "far", TypeError),
        ("operation", "mirror", "flip", ValueError),
        ("operation", "mirror", 123, TypeError),
        ("rotation", 3, "three", TypeError),
        ("rotation", 3, 1, ValueError),
    ],
)
def test_symmetrization_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = SymmetrizationConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("distance", 3.0),
        ("operation", "mirror"),
        ("rotation", None),
    ],
)
def test_symmetrization_config_default_values(property_name, initial_value):
    """
    Test default values for SymmetrizationConfig properties.
    """
    config = SymmetrizationConfig()
    assert getattr(config, property_name) == initial_value
