"""
Contains the configuration Class for the program.
"""

from __future__ import annotations

from pathlib import Path
from abc import ABC, abstractmethod
import toml

from ..molecules import PSE_NUMBERS


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

    def __init__(self: GeneralConfig) -> None:
        self._verbosity: int = 1
        self._max_cycles: int = 100
        self._print_config: bool = False
        self._parallel: int = 1
        self._num_molecules: int = 1
        self._postprocess: bool = False
        self._write_xyz: bool = True

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
        if verbosity not in [-1, 0, 1, 2, 3]:
            raise ValueError("Verbosity can only be -1, 0, 1, 2, or 3.")
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

    @property
    def num_molecules(self):
        """
        Get the number of molecules.
        """
        return self._num_molecules

    @num_molecules.setter
    def num_molecules(self, num_molecules: int):
        """
        Set the number of molecules.
        """
        if not isinstance(num_molecules, int):
            raise TypeError("Number of molecules should be an integer.")
        if num_molecules < 1:
            raise ValueError("Number of molecules should be greater than 0.")
        self._num_molecules = num_molecules

    @property
    def postprocess(self):
        """
        Get the postprocess flag.
        """
        return self._postprocess

    @postprocess.setter
    def postprocess(self, postprocess: bool):
        """
        Set the postprocess flag.
        """
        if not isinstance(postprocess, bool):
            raise TypeError("Postprocess should be a boolean.")
        self._postprocess = postprocess

    @property
    def write_xyz(self):
        """
        Get the write xyz flag.
        """
        return self._write_xyz

    @write_xyz.setter
    def write_xyz(self, write_xyz: bool):
        """
        Set the write xyz flag.
        """
        if not isinstance(write_xyz, bool):
            raise TypeError("Write xyz should be a boolean.")
        self._write_xyz = write_xyz


class GenerateConfig(BaseConfig):
    """
    Configuration for the "generate" section responsible for setting up an initial Molecule type.
    """

    def __init__(self: GenerateConfig) -> None:
        self._min_num_atoms: int = 2
        self._max_num_atoms: int = 100
        self._init_coord_scaling: float = 3.0
        self._dist_threshold: float = 1.2
        self._increase_scaling_factor: float = 1.3
        self._element_composition: dict[int, tuple[int | None, int | None]] = {}
        self._forbidden_elements: list[int] | None = None

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

    @property
    def init_coord_scaling(self):
        """
        Get the initial coordinate scaling.
        """
        return self._init_coord_scaling

    @init_coord_scaling.setter
    def init_coord_scaling(self, init_coord_scaling: float):
        """
        Set the initial coordinate scaling.
        """
        if not isinstance(init_coord_scaling, float):
            raise TypeError("Initial coordinate scaling should be a float.")
        if init_coord_scaling <= 0:
            raise ValueError("Initial coordinate scaling should be greater than 0.")
        self._init_coord_scaling = init_coord_scaling

    @property
    def dist_threshold(self):
        """
        Get the distance threshold.
        """
        return self._dist_threshold

    @dist_threshold.setter
    def dist_threshold(self, dist_threshold: float):
        """
        Set the distance threshold.
        """
        if not isinstance(dist_threshold, float):
            raise TypeError("Distance threshold should be a float.")
        if dist_threshold <= 0:
            raise ValueError("Distance threshold should be greater than 0.")
        self._dist_threshold = dist_threshold

    @property
    def increase_scaling_factor(self):
        """
        Get the increase scaling factor.
        """
        return self._increase_scaling_factor

    @increase_scaling_factor.setter
    def increase_scaling_factor(self, increase_scaling_factor: float):
        """
        Set the increase scaling factor.
        """
        if not isinstance(increase_scaling_factor, float):
            raise TypeError("Increase scaling factor should be a float.")
        if increase_scaling_factor <= 1:
            raise ValueError("Increase scaling factor should be greater than 1.")
        self._increase_scaling_factor = increase_scaling_factor

    @property
    def element_composition(self):
        """
        Return the element composition.
        """
        return self._element_composition

    @element_composition.setter
    def element_composition(self, composition_str):
        """
        Parses the element_composition string and stores the parsed data
        in the _element_composition dictionary.

        Format: "C:2-10, H:10-20, O:1-5, N:1-*"
        """

        if not isinstance(composition_str, str):
            raise TypeError("Element composition should be a string.")
        if not composition_str:
            return

        element_dict = {}
        elements = composition_str.split(",")
        # remove leading and trailing whitespaces
        elements = [element.strip() for element in elements]

        for element in elements:
            element_type, range_str = element.split(":")
            min_count, max_count = range_str.split("-")
            element_number = PSE_NUMBERS.get(element_type.lower(), None) - 1

            # Convert counts, handle wildcard '*'
            min_count = None if min_count == "*" else int(min_count)
            max_count = None if max_count == "*" else int(max_count)
            if (
                min_count is not None
                and max_count is not None
                and min_count > max_count
            ):
                raise ValueError(
                    f"Minimum count ({min_count}) is larger than maximum count ({max_count})."
                )

            element_dict[element_number] = (min_count, max_count)

        self._element_composition = element_dict

    @property
    def forbidden_elements(self):
        """
        Return sorted list of forbidden elements.
        """
        return self._forbidden_elements

    @forbidden_elements.setter
    def forbidden_elements(self: GenerateConfig, forbidden_str: str) -> None:
        """
        Parses the forbidden_elements string and stores the parsed data
        in the _forbidden_elements set.

        Format: "57-71, 8, 1" or "19-*"
        """
        # if string is empty or None, set to None
        if not forbidden_str:
            self._forbidden_elements = None
            return
        forbidden_set: set[int] = set()
        elements = forbidden_str.split(",")
        elements = [element.strip() for element in elements]

        for item in elements:
            if "-" in item:
                start: str | int
                end: str | int
                start, end = item.split("-")
                if end == "*" and start == "*":
                    raise ValueError("Both start and end cannot be wildcard '*'.")
                if end == "*":
                    end = 103  # Set to a the maximum atomic number in mindlessgen
                if start == "*":
                    start = 0
                forbidden_set.update(
                    range(int(start) - 1, int(end))
                )  # Subtract 1 to convert to 0-based indexing
            else:
                forbidden_set.add(
                    int(item) - 1
                )  # Subtract 1 to convert to 0-based indexing

        self._forbidden_elements = sorted(list(forbidden_set))


class RefineConfig(BaseConfig):
    """
    Configuration class for refinement settings.
    """

    def __init__(self: RefineConfig) -> None:
        self._max_frag_cycles: int = 100
        self._engine: str = "xtb"
        self._hlgap: float = 0.5
        self._debug: bool = False

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
            raise TypeError("Refinement engine should be a string.")
        if engine not in ["xtb", "orca"]:
            raise ValueError("Refinement engine can only be xtb or orca.")
        self._engine = engine

    @property
    def hlgap(self):
        """
        Get the minimum HOM
        """
        return self._hlgap

    @hlgap.setter
    def hlgap(self, hlgap: float):
        """
        Set the minimum HOM
        """
        if not isinstance(hlgap, float):
            raise TypeError("Minimum HL gap should be a float.")
        if hlgap < 0:
            raise ValueError("Minimum HL gap should be greater than 0.")
        self._hlgap = hlgap

    @property
    def debug(self):
        """
        Get the debug flag for refinement.
        """
        return self._debug

    @debug.setter
    def debug(self, debug: bool):
        """
        Set the debug flag for refinement.
        """
        if not isinstance(debug, bool):
            raise TypeError("Debug should be a boolean.")
        self._debug = debug


class PostProcessConfig(BaseConfig):
    """
    Configuration class for post-processing settings.
    """

    def __init__(self: PostProcessConfig) -> None:
        self._engine: str = "orca"
        self._opt_cycles: int | None = None
        self._optimize: bool = False
        self._debug: bool = False

    def get_identifier(self) -> str:
        return "postprocess"

    @property
    def engine(self):
        """
        Get the postprocess engine.
        """
        return self._engine

    @engine.setter
    def engine(self, engine: str):
        """
        Set the postprocess engine.
        """
        if not isinstance(engine, str):
            raise TypeError("Postprocess engine should be a string.")
        if engine not in ["xtb", "orca"]:
            raise ValueError("Postprocess engine can only be xtb or orca.")
        self._engine = engine

    @property
    def optimize(self):
        """
        Get the optimization flag for post-processing.
        """
        return self._optimize

    @optimize.setter
    def optimize(self, optimize: bool):
        """
        Set the optimization flag for post-processing.
        """
        if not isinstance(optimize, bool):
            raise TypeError("Optimize should be a boolean.")
        self._optimize = optimize

    @property
    def opt_cycles(self):
        """
        Get the optimization cycles for post-processing.
        """
        return self._opt_cycles

    @opt_cycles.setter
    def opt_cycles(self, opt_cycles: int):
        """
        Set the optimization cycles for post-processing.
        """
        if not isinstance(opt_cycles, int):
            raise TypeError("Optimization cycles should be an integer.")
        if opt_cycles < 0:
            raise ValueError("Optimization cycles should be 0 or greater.")
        self._opt_cycles = opt_cycles

    @property
    def debug(self):
        """
        Get the debug flag for post-processing.
        """
        return self._debug

    @debug.setter
    def debug(self, debug: bool):
        """
        Set the debug flag for post-processing.
        """
        if not isinstance(debug, bool):
            raise TypeError("Debug should be a boolean.")
        self._debug = debug


class XTBConfig(BaseConfig):
    """
    Configuration class for XTB.
    """

    def __init__(self: XTBConfig) -> None:
        self._xtb_path: str | Path = "xtb"
        self._level: int = 2

    def get_identifier(self) -> str:
        return "xtb"

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

    @property
    def level(self):
        """
        Get the GFN<n>-xTB level.
        """
        return self._level

    @level.setter
    def level(self, level: int):
        """
        Set the GFN<n>-xTB level.
        """
        if not isinstance(level, int):
            raise TypeError("xtb level should be an integer.")
        self._level = level


class ORCAConfig(BaseConfig):
    """
    Configuration class for ORCA.
    """

    def __init__(self: ORCAConfig) -> None:
        self._orca_path: str | Path = "orca"
        self._functional: str = "PBE"
        self._basis: str = ""
        self._gridsize: int = 1
        self._scf_cycles: int = 100

    def get_identifier(self) -> str:
        return "orca"

    @property
    def orca_path(self):
        """
        Get the orca path.
        """
        return self._orca_path

    @orca_path.setter
    def orca_path(self, orca_path: str | Path):
        """
        Set the orca path.
        """
        if not isinstance(orca_path, str | Path):
            raise TypeError("orca_path should be a string or Path.")
        self._orca_path = orca_path

    @property
    def functional(self):
        """
        Get the ORCA functional/method.
        """
        return self._functional

    @functional.setter
    def functional(self, functional: str):
        """
        Set the ORCA functional/method.
        """
        if not isinstance(functional, str):
            raise TypeError("Functional should be a string.")
        self._functional = functional

    @property
    def basis(self):
        """
        Get the ORCA basis set.
        """
        return self._basis

    @basis.setter
    def basis(self, basis: str):
        """
        Set the ORCA basis set.
        """
        if not isinstance(basis, str):
            raise TypeError("Basis should be a string.")
        self._basis = basis

    @property
    def gridsize(self):
        """
        Get the ORCA gridsize for numerical integration.
        """
        return self._gridsize

    @gridsize.setter
    def gridsize(self, gridsize: int):
        """
        Set the ORCA gridsize for numerical integration.
        """
        if not isinstance(gridsize, int):
            raise TypeError("Gridsize should be an integer.")
        if gridsize not in [1, 2, 3]:
            raise ValueError(
                "Gridsize should be 1, 2, or 3. (-> DEFGRID1, DEFGRID2, DEFGRID3)"
            )
        self._gridsize = gridsize

    @property
    def scf_cycles(self):
        """
        Get the maximum number of SCF cycles.
        """
        return self._scf_cycles

    @scf_cycles.setter
    def scf_cycles(self, max_scf_cycles: int):
        """
        Set the maximum number of SCF cycles.
        """
        if not isinstance(max_scf_cycles, int):
            raise TypeError("Max SCF cycles should be an integer.")
        if max_scf_cycles < 1:
            raise ValueError("Max SCF cycles should be greater than 0.")
        self._scf_cycles = max_scf_cycles


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
        self.postprocess = PostProcessConfig()
        self.generate = GenerateConfig()

        if config_file:
            self.load_from_toml(config_file)

    def get_all_identifiers(self):
        """
        Returns the identifiers of all subconfiguration classes, e.g. "orca", "refinement", ...
        """
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
        all_identifiers = self.get_all_identifiers()
        for key in config_dict:
            if key not in all_identifiers:
                raise KeyError(f"Unknown key in configuration file: {key}")

        for sub_config in all_identifiers:
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
