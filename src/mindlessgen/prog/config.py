"""
Contains the configuration Class for the program.
"""

from pathlib import Path
import toml
from abc import ABC, abstractmethod


# abstract base class for configuration
class BaseConfig(ABC):
    """
    Abstract base class for configuration settings.
    """

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_identifier(self) -> str:
        """
        Get the identifier of the configuration.
        """


class GeneralConfig(BaseConfig):
    """
    Configuration class for general settings.
    """

    def __init__(self):
        self._verbosity: int = 1
        self._max_cycles: int = 100
        self._engine: str = "xtb"
        self._print_config: bool = False
        self._parallel: int = 1

    def get_identifier(self) -> str:
        return "general"

    @property
    def verbosity(self):
        """
        Get the verbosity level.
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        """
        Set the verbosity level.
        """
        if not isinstance(verbosity, int):
            raise TypeError("Verbosity should be an integer.")
        if verbosity not in [0, 1, 2, 3]:
            raise ValueError("Verbosity can only be 0, 1, 2, or 3.")
        self._verbosity = verbosity

    @property
    def max_cycles(self):
        """
        Get the maximum number of cycles.
        """
        return self._max_cycles

    @max_cycles.setter
    def max_cycles(self, max_cycles: int):
        """
        Set the maximum number of cycles.
        """
        if not isinstance(max_cycles, int):
            raise TypeError("Max cycles should be an integer.")
        if max_cycles < 1:
            raise ValueError("Max cycles should be greater than 0.")
        self._max_cycles = max_cycles

    @property
    def engine(self):
        """
        Get the engine.
        """
        return self._engine

    @engine.setter
    def engine(self, engine: str):
        """
        Set the engine.
        """
        if not isinstance(engine, str):
            raise TypeError("Engine should be a string.")
        if engine not in ["xtb", "orca"]:
            raise ValueError("Engine can only be xtb or orca.")
        self._engine = engine

    @property
    def print_config(self):
        """
        Get the print config flag.
        """
        return self._print_config

    @print_config.setter
    def print_config(self, print_config: bool):
        """
        Set the print config flag.
        """
        if not isinstance(print_config, bool):
            raise TypeError("Print config should be a boolean.")
        self._print_config = print_config

    @property
    def parallel(self):
        """
        Get the parallel flag.
        """
        return self._parallel

    @parallel.setter
    def parallel(self, parallel: int):
        """
        Set the parallel flag.
        """
        if not isinstance(parallel, int):
            raise TypeError("Parallel should be an integer.")
        if parallel < 1:
            raise ValueError("Parallel should be greater than 0.")
        self._parallel = parallel


class GenerateConfig(BaseConfig):
    def __init__(self):
        self._min_num_atoms: int = 2
        self._max_num_atoms: int = 100

    def get_identifier(self) -> str:
        return "generate"

    @property
    def min_num_atoms(self):
        """
        Get the minimum number of atoms.
        """
        return self._min_num_atoms

    @min_num_atoms.setter
    def min_num_atoms(self, min_num_atoms: int):
        """
        Set the minimum number of atoms.
        """
        if not isinstance(min_num_atoms, int):
            raise TypeError("Min num atoms should be an integer.")
        if min_num_atoms < 1:
            raise ValueError("Min num atoms should be greater than 0.")
        self._min_num_atoms = min_num_atoms

    @property
    def max_num_atoms(self):
        """
        Get the maximum number of atoms.
        """
        return self._max_num_atoms

    @max_num_atoms.setter
    def max_num_atoms(self, max_num_atoms: int):
        """
        Set the maximum number of atoms.
        """
        if not isinstance(max_num_atoms, int):
            raise TypeError("Max num atoms should be an integer.")
        if max_num_atoms < 1:
            raise ValueError("Max num atoms should be greater than 0.")
        self._max_num_atoms = max_num_atoms


class RefineConfig(BaseConfig):
    """
    Configuration class for refinement settings.
    """

    def __init__(self):
        self._max_frag_cycles: int = 100

    def get_identifier(self) -> str:
        return "refine"

    @property
    def max_frag_cycles(self):
        """
        Get the maximum number of fragment cycles.
        """
        return self._max_frag_cycles

    @max_frag_cycles.setter
    def max_frag_cycles(self, max_frag_cycles: int):
        """
        Set the maximum number of fragment cycles.
        """
        if not isinstance(max_frag_cycles, int):
            raise TypeError("Max fragment cycles should be an integer.")
        if max_frag_cycles < 1:
            raise ValueError("Max fragment cycles should be greater than 0.")
        self._max_frag_cycles = max_frag_cycles


class XTBConfig(BaseConfig):
    """
    Configuration class for XTB.
    """

    def __init__(self):
        self._xtb_option: str = "dummy"
        self._xtb_path: str | Path = "xtb"

    def get_identifier(self) -> str:
        return "xtb"

    @property
    def xtb_option(self):
        """
        Get the xtb option.
        """
        return self._xtb_option

    @xtb_option.setter
    def xtb_option(self, xtb_option: str):
        """
        Set the xtb option.
        """
        if not isinstance(xtb_option, str):
            raise TypeError("xtb_option should be a string.")
        self._xtb_option = xtb_option

    @property
    def xtb_path(self):
        """
        Get the xtb path.
        """
        return self._xtb_path

    @xtb_path.setter
    def xtb_path(self, xtb_path: str | Path):
        """
        Set the xtb path.
        """
        if not isinstance(xtb_path, str | Path):
            raise TypeError("xtb_path should be a string.")
        self._xtb_path = xtb_path


class ORCAConfig(BaseConfig):
    """
    Configuration class for ORCA.
    """

    def __init__(self):
        self._orca_option: str = "dummy"

    def get_identifier(self) -> str:
        return "orca"

    @property
    def orca_option(self):
        """
        Get the orca option.
        """
        return self._orca_option

    @orca_option.setter
    def orca_option(self, orca_option: str):
        """
        Set the orca option.
        """
        if not isinstance(orca_option, str):
            raise TypeError("orca_option should be a string.")
        self._orca_option = orca_option


class ConfigManager:
    """
    Overall configuration manager for the program.
    """

    def __init__(self, config_file: str | Path | None = None):
        """
        Initialize configuration sections with default values
        """
        self.general = GeneralConfig()
        self.xtb = XTBConfig()
        self.orca = ORCAConfig()
        self.refine = RefineConfig()
        self.generate = GenerateConfig()

        if config_file:
            self.load_from_toml(config_file)

    def get_all_identifiers(self):
        identifiers = []
        for attr_name in dir(self):
            attr_value = getattr(self, attr_name)
            # Check if the attribute is an instance of BaseConfig
            if isinstance(attr_value, BaseConfig):
                identifiers.append(attr_value.get_identifier())
        return identifiers

    def load_from_toml(self, config_file: str | Path) -> None:
        """
        Load configuration from TOML file that is structured as follows:
        [general]
        verbosity = 1
        max_cycles = 100
        engine = "xtb"

        [xtb]
        xtb_option = "opt"

        [orca]
        orca_option = "opt"

        Arguments:
            config_file (str): Path to the configuration file

        """
        # Load the configuration file
        config_data = toml.load(config_file)
        self.load_from_dict(config_data)

    def load_from_dict(self, config_dict: dict) -> None:
        """
        Load configuration from a dictionary structured as follows:
        {
            "general": {
                "verbosity": 1,
                "max_cycles": 100,
                "input": "input.xyz"
            },
            "xtb": {
                "xtb_option": "opt"
            },
            "orca": {
                "orca_option": "opt"
            }
        }

        Arguments:
            config_dict (dict): Dictionary containing the configuration
        """
        # Check for unknown keys
        for key in config_dict:
            if key not in self.get_all_identifiers():
                raise KeyError(f"Unknown key in configuration file: {key}")

        for sub_config in self.get_all_identifiers():
            if sub_config not in config_dict:
                continue
            for config_key, config_value in config_dict[sub_config].items():
                # check if config_value is not None and if the attribute exists
                if config_value is not None and hasattr(
                    getattr(self, sub_config), config_key
                ):
                    setattr(getattr(self, sub_config), config_key, config_value)

    def __str__(self) -> str:
        """
        Automated method to display the current configuration.
        """
        configstr = ""
        for attr_name in dir(self):
            attr_value = getattr(self, attr_name)
            if isinstance(attr_value, BaseConfig):
                configstr += (
                    f"{attr_value.get_identifier().capitalize()} configuration:\n"
                )
                for key, value in attr_value.__dict__.items():
                    configstr += (
                        f"{key[1:]:>30}:   {value}\n"  # Skip the leading underscore
                    )
                configstr += "\n"
        return configstr
