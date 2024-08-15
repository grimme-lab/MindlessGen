import pytest
from mindlessgen.prog import GeneralConfig  # type: ignore


@pytest.mark.parametrize(
    "property_name, valid_value, invalid_value, expected_exception",
    [
        ("verbosity", 2, "high", TypeError),
        ("verbosity", 2, 4, ValueError),
        ("max_cycles", 50, -1, ValueError),
        ("max_cycles", 50, "100", TypeError),
        ("engine", "orca", "mopac", ValueError),
        ("engine", "orca", 123, TypeError),
        ("min_num_atoms", 5, 0, ValueError),
        ("min_num_atoms", 5, "two", TypeError),
        ("max_num_atoms", 80, -10, ValueError),
        ("max_num_atoms", 80, None, TypeError),
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


@pytest.mark.parametrize(
    "property_name, initial_value",
    [
        ("verbosity", 1),
        ("max_cycles", 100),
        ("engine", "xtb"),
        ("min_num_atoms", 2),
        ("max_num_atoms", 100),
        ("print_config", False),
        ("parallel", 1),
    ],
)
def test_general_config_default_values(property_name, initial_value):
    config = GeneralConfig()
    assert getattr(config, property_name) == initial_value
