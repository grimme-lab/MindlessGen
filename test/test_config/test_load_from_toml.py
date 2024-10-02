import pytest
from pathlib import Path
from mindlessgen.prog import ConfigManager  # type: ignore


@pytest.fixture
def toml_file_path():
    # Replace this with the path to your actual TOML file
    # Specify the path to the actual TOML file
    testsdir = Path(__file__).resolve().parents[1]
    toml_file = testsdir / "fixtures/example_config.toml"
    return Path(toml_file).resolve()


@pytest.fixture
def config_manager(toml_file_path):
    # Ensure the TOML file exists before proceeding
    assert toml_file_path.exists(), f"TOML file not found at: {toml_file_path}"

    # Load the configuration from the TOML file
    return ConfigManager(config_file=toml_file_path)


def test_load_general_config(config_manager):
    # Verify 'general' settings
    assert config_manager.general.verbosity == 1
    assert config_manager.general.parallel == 1
    assert config_manager.general.max_cycles == 100
    assert config_manager.general.num_molecules == 1
    assert config_manager.general.postprocess is False


def test_load_generate_config(config_manager):
    # Verify 'generate' settings
    assert config_manager.generate.min_num_atoms == 2
    assert config_manager.generate.max_num_atoms == 100
    assert config_manager.generate.init_coord_scaling == 3.0
    assert config_manager.generate.increase_scaling_factor == 1.3
    assert config_manager.generate.element_composition == {
        5: (2, 10),  # Carbon (C)
        0: (10, 20),  # Hydrogen (H)
        7: (1, 5),  # Oxygen (O)
        6: (1, None),  # Nitrogen (N)
    }
    assert config_manager.generate.forbidden_elements == list(range(56, 71))


def test_load_refine_config(config_manager):
    # Verify 'refine' settings
    assert config_manager.refine.max_frag_cycles == 100
    assert config_manager.refine.engine == "xtb"


def test_load_postprocess_config(config_manager):
    # Verify 'postprocess' settings
    assert config_manager.postprocess.engine == "xtb"


def test_load_xtb_config(config_manager):
    # Verify 'xtb' settings
    assert config_manager.xtb.xtb_path == "/path/to/xtb"


def test_load_orca_config(config_manager):
    # Verify 'orca' settings
    assert config_manager.orca.orca_path == "/path/to/orca"
