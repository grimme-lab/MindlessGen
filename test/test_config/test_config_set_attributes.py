import pytest
from mindlessgen.prog import (  # type: ignore
    GeneralConfig,
    GenerateConfig,
    RefineConfig,
    XTBConfig,
)


@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("verbosity", 2, "high", TypeError),
        ("verbosity", 2, 4, ValueError),
        ("max_cycles", 50, -1, ValueError),
        ("max_cycles", 50, "100", TypeError),
        ("engine", "orca", "mopac", ValueError),
        ("engine", "orca", 123, TypeError),
        ("print_config", True, "yes", TypeError),
        ("parallel", 4, 0, ValueError),
        ("parallel", 4, "four", TypeError),
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


# create a similar test for GenerateConfig
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


# create a similar test for RefineConfig
@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("max_frag_cycles", 100, -1, ValueError),
        ("max_frag_cycles", 100, "100", TypeError),
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
        ("verbosity", 1),
        ("max_cycles", 100),
        ("engine", "xtb"),
        ("print_config", False),
        ("parallel", 1),
    ],
)
def test_general_config_default_values(property_name, initial_value):
    config = GeneralConfig()
    assert getattr(config, property_name) == initial_value


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("min_num_atoms", 2),
        ("max_num_atoms", 100),
        ("init_coord_scaling", 3.0),
        ("dist_threshold", 1.2),
        ("increase_scaling_factor", 1.3),
    ],
)
def test_generate_config_default_values(property_name, initial_value):
    config = GenerateConfig()
    assert getattr(config, property_name) == initial_value


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("max_frag_cycles", 100),
    ],
)
def test_refine_config_default_values(property_name, initial_value):
    config = RefineConfig()
    assert getattr(config, property_name) == initial_value


# Generate tests for XTBConfig
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
