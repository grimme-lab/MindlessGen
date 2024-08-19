import pytest
from mindlessgen.prog import (  # type: ignore
    GeneralConfig,
    GenerateConfig,
    RefineConfig,
    XTBConfig,
    ORCAConfig,
    PostProcessConfig,
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
        ("max_cycles", 100),
        ("print_config", False),
        ("parallel", 1),
        ("num_molecules", 1),
        ("postprocess", False),
    ],
)
def test_general_config_default_values(property_name, initial_value):
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
        ("dist_threshold", 1.5, -1.0, ValueError),
        ("dist_threshold", 1.5, "1.5", TypeError),
        ("increase_scaling_factor", 1.1, 0.0, ValueError),
        ("increase_scaling_factor", 1.1, "1.1", TypeError),
    ],
)
def test_generate_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
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
        ("forbidden_elements", "6,1", [0, 5]),
        (
            "forbidden_elements",
            "86-*",
            [85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102],
        ),
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
        ("min_num_atoms", 2),
        ("max_num_atoms", 100),
        ("init_coord_scaling", 3.0),
        ("dist_threshold", 1.2),
        ("increase_scaling_factor", 1.3),
        ("element_composition", {}),
        ("forbidden_elements", None),
    ],
)
def test_generate_config_default_values(property_name, initial_value):
    config = GenerateConfig()
    assert getattr(config, property_name) == initial_value


# Tests for RefineConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("max_frag_cycles", 100, -1, ValueError),
        ("max_frag_cycles", 100, "100", TypeError),
        ("engine", "xtb", 123, TypeError),
        ("engine", "xtb", "g16", ValueError),
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
        ("max_frag_cycles", 100),
        ("engine", "xtb"),
    ],
)
def test_refine_config_default_values(property_name, initial_value):
    config = RefineConfig()
    assert getattr(config, property_name) == initial_value


# Tests for XTBConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("xtb_path", "path/to/xtb", 123, TypeError),
        ("xtb_path", "path/to/xtb", None, TypeError),
    ],
)
def test_xtb_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
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
    ],
)
def test_xtb_config_default_values(property_name, initial_value):
    config = XTBConfig()
    assert getattr(config, property_name) == initial_value


# Tests for ORCAConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("orca_path", "path/to/orca", 123, TypeError),
        ("orca_path", "path/to/orca", None, TypeError),
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
    ],
)
def test_orca_config_default_values(property_name, initial_value):
    config = ORCAConfig()
    assert getattr(config, property_name) == initial_value


# Tests for PostProcessConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("engine", "orca", 123, TypeError),
        ("engine", "orca", "g16", ValueError),
    ],
)
def test_postprocess_config_property_setters(
    property_name, valid_value, invalid_value, expected_exception
):
    config = PostProcessConfig()

    # Test valid value
    setattr(config, property_name, valid_value)
    assert getattr(config, property_name) == valid_value

    # Test invalid value
    with pytest.raises(expected_exception):
        setattr(config, property_name, invalid_value)


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("engine", "orca"),
    ],
)
def test_postprocess_config_default_values(property_name, initial_value):
    config = PostProcessConfig()
    assert getattr(config, property_name) == initial_value
