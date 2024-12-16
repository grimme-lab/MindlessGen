"""
Contains the configuration Class for the program.
"""

from __future__ import annotations

from pathlib import Path
from abc import ABC, abstractmethod
import warnings
import multiprocessing as mp

import numpy as np
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
        self._max_cycles: int = 200
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

    def check_config(self, verbosity: int = 1) -> None:
        ### GeneralConfig checks ###
        # lower number of the available cores and the configured parallelism
        num_cores = min(mp.cpu_count(), self.parallel)
        if self.parallel > mp.cpu_count() and verbosity > -1:
            warnings.warn(
                f"Number of cores requested ({self.parallel}) is greater "
                + f"than the number of available cores ({mp.cpu_count()})."
                + f"Using {num_cores} cores instead."
            )

        if num_cores > 1 and verbosity > 0:
            # raise warning that parallelization will disable verbosity
            warnings.warn(
                "Parallelization will disable verbosity during iterative search. "
                + "Set '--verbosity 0' or '-P 1' to avoid this warning, or simply ignore it."
            )


class GenerateConfig(BaseConfig):
    """
    Configuration for the "generate" section responsible for setting up an initial Molecule type.
    """

    def __init__(self: GenerateConfig) -> None:
        self._min_num_atoms: int = 5
        self._max_num_atoms: int = 10
        self._init_coord_scaling: float = 3.0
        self._increase_scaling_factor: float = 1.1
        self._element_composition: dict[int, tuple[int | None, int | None]] = {}
        self._forbidden_elements: list[int] | None = None
        self._scale_fragment_detection: float = 1.25
        self._scale_minimal_distance: float = 0.8
        self._contract_coords: bool = True
        self._molecular_charge: int | None = None
        self._fixed_composition: bool = False

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
    def element_composition(
        self, composition: None | str | dict[int, tuple[int | None, int | None]]
    ) -> None:
        """
        If composition_str: str, it should be a string with the format:
            Parses the element_composition string and stores the parsed data
            in the _element_composition dictionary.
            Format: "C:2-10, H:10-20, O:1-5, N:1-*"
        If composition_str: dict, it should be a dictionary with integer keys and tuple values. Will be stored as is.

        Arguments:
            composition_str (str): String with the element composition
            composition_str (dict): Dictionary with integer keys and tuple values
        Raises:
            TypeError: If composition_str is not a string or a dictionary
            AttributeError: If the element is not found in the periodic table
            ValueError: If the minimum count is larger than the maximum count
        Returns:
            None
        """

        if not composition:
            return
        if isinstance(composition, dict):
            for key, value in composition.items():
                if (
                    not isinstance(key, int)
                    or not isinstance(value, tuple)
                    or len(value) != 2
                    or not all(isinstance(val, int) or val is None for val in value)
                ):
                    raise TypeError(
                        "Element composition dictionary should be a dictionary with integer keys and tuple values (int, int)."
                    )
            self._element_composition = composition
            return
        if not isinstance(composition, str):
            raise TypeError(
                "Element composition should be a string (will be parsed) or "
                + "a dictionary with integer keys and tuple values."
            )

        element_dict: dict[int, tuple[int | None, int | None]] = {}
        elements = composition.split(",")
        # remove leading and trailing whitespaces
        elements = [element.strip() for element in elements]

        min_count: int | str | None
        max_count: int | str | None
        for element in elements:
            element_type, range_str = element.split(":")
            min_count, max_count = range_str.split("-")
            element_number = PSE_NUMBERS.get(element_type.lower(), None)
            if element_number is None:
                raise AttributeError(
                    f"Element {element_type} not found in the periodic table."
                )
            # correct for 1- vs. 0-based indexing
            element_number = element_number - 1

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
    def forbidden_elements(
        self: GenerateConfig, forbidden: None | str | list[int]
    ) -> None:
        """
        If forbidden: str:
            Parses the forbidden_elements string and stores the parsed data
            in the _forbidden_elements set.
            Format: "57-71, 8, 1" or "19-*"
        If forbidden: list:
            Stores the forbidden elements as is.

        Arguments:
            forbidden (str): String with the forbidden elements
            forbidden (list): List with integer values
        Raises:
            TypeError: If forbidden is not a string or a list of integers
            ValueError: If both start and end are wildcard '*'
        Returns:
            None
        """
        # if string is empty or None, set to None
        if not forbidden:
            self._forbidden_elements = None
            return
        if isinstance(forbidden, list):
            if all(isinstance(elem, int) for elem in forbidden):
                self._forbidden_elements = sorted(forbidden)
                return
            raise TypeError("Forbidden elements should be a list of integers.")
        if not isinstance(forbidden, str):
            raise TypeError(
                "Forbidden elements should be a string or a list of integers."
            )
        forbidden_set: set[int] = set()
        elements = forbidden.split(",")
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

    @property
    def scale_fragment_detection(self):
        """
        Get the scaling factor for the fracment detection based on the van der Waals radii.
        """
        return self._scale_fragment_detection

    @scale_fragment_detection.setter
    def scale_fragment_detection(self, scale_fragment_detection: float):
        """
        Set the scaling factor for van der Waals radii.
        """
        if not isinstance(scale_fragment_detection, float):
            raise TypeError("Scale van der Waals radii should be a float.")
        if scale_fragment_detection <= 0:
            raise ValueError("Scale van der Waals radii should be greater than 0.")
        self._scale_fragment_detection = scale_fragment_detection

    @property
    def scale_minimal_distance(self):
        """
        Get the scaling factor for minimal distance between two atoms.
        """
        return self._scale_minimal_distance

    @scale_minimal_distance.setter
    def scale_minimal_distance(self, scale_minimal_distance: float):
        """
        Set the scaling factor for minimal distance between two atoms.
        """
        if not isinstance(scale_minimal_distance, float):
            raise TypeError("Scale minimal distance should be a float.")
        if scale_minimal_distance <= 0:
            raise ValueError("Scale minimal distance should be greater than 0.")
        self._scale_minimal_distance = scale_minimal_distance

    @property
    def contract_coords(self):
        """
        Get the contract_coords flag.
        """
        return self._contract_coords

    @contract_coords.setter
    def contract_coords(self, contract_coords: bool):
        """
        Set the contract_coords flag.
        """
        if not isinstance(contract_coords, bool):
            raise TypeError("Contract coords should be a boolean.")
        self._contract_coords = contract_coords

    @property
    def molecular_charge(self):
        """
        Get the molecular_charge.
        """
        return self._molecular_charge

    @molecular_charge.setter
    def molecular_charge(self, molecular_charge: str | int):
        """
        Set the molecular_charge.
        """
        if isinstance(molecular_charge, str):
            if molecular_charge.lower() == "none" or molecular_charge == "":
                self._molecular_charge = None
            else:
                self._molecular_charge = int(molecular_charge)
        elif isinstance(molecular_charge, int):
            self._molecular_charge = molecular_charge
        else:
            raise TypeError("Molecular charge should be a string or an integer.")

    @property
    def fixed_composition(self):
        """
        Get the fixed_composition flag.
        """
        return self._fixed_composition

    @fixed_composition.setter
    def fixed_composition(self, fixed_composition: bool):
        """
        Set the fixed_composition flag.
        """
        if not isinstance(fixed_composition, bool):
            raise TypeError("Fixed composition should be a boolean.")
        self._fixed_composition = fixed_composition

    def check_config(self, verbosity: int = 1) -> None:
        """
        GenerateConfig checks for any incompatibilities that are imaginable
        """
        # - Check if the minimum number of atoms is smaller than the maximum number of atoms
        if self.min_num_atoms is not None and self.max_num_atoms is not None:
            if self.min_num_atoms > self.max_num_atoms:
                raise ValueError(
                    "The minimum number of atoms is larger than the maximum number of atoms."
                )
        # - Check if the summed number of minimally required atoms from cfg.element_composition
        #   is larger than the maximum number of atoms
        if self.max_num_atoms is not None:
            if (
                np.sum(
                    [
                        self.element_composition.get(i, (0, 0))[0]
                        if self.element_composition.get(i, (0, 0))[0] is not None
                        else 0
                        for i in self.element_composition
                    ]
                )
                > self.max_num_atoms
            ):
                raise ValueError(
                    "The summed number of minimally required atoms "
                    + "from the fixed composition is larger than the maximum number of atoms."
                )
        if self.fixed_composition:
            # - Check if all defintions in cfg.element_composition
            #   are completely fixed (min and max are equal)
            for elem, count_range in self.element_composition.items():
                # check if for all entries: min and max are both not None, and if min and max are equal.
                if (
                    (count_range[0] is None)
                    or (count_range[1] is None)
                    or (count_range[0] != count_range[1])
                ):
                    raise ValueError(
                        f"Element {elem} is not completely fixed in the element composition. "
                        + "Usage together with fixed_composition is not possible."
                    )
            # Check if the summed number of fixed atoms
            # is within the defined overall limits
            sum_fixed_atoms = np.sum(
                [
                    self.element_composition.get(i, (0, 0))[0]
                    for i in self.element_composition
                ]
            )
            if self.min_num_atoms is not None and not (
                self.min_num_atoms <= sum_fixed_atoms
            ):
                raise ValueError(
                    "The summed number of fixed atoms from the fixed composition "
                    + "is not within the range of min_num_atoms and max_num_atoms."
                )


class RefineConfig(BaseConfig):
    """
    Configuration class for refinement settings.
    """

    def __init__(self: RefineConfig) -> None:
        self._max_frag_cycles: int = 10
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
        self._opt_cycles: int | None = 5
        self._optimize: bool = True
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
        if engine not in ["xtb", "orca", "gp3"]:
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
        self._basis: str = "def2-SVP"
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

    def check_config(self, verbosity: int = 1) -> None:
        """
        Checks ConfigClass for any incompatibilities that are imaginable
        """

        ### Config-specific checks ###
        for attr_name in dir(self):
            attr_value = getattr(self, attr_name)
            if isinstance(attr_value, BaseConfig):
                if hasattr(attr_value, "check_config"):
                    attr_value.check_config(verbosity)

        ### Overlapping checks ###
        num_cores = min(mp.cpu_count(), self.general.parallel)
        if num_cores > 1 and self.postprocess.debug and verbosity > -1:
            # raise warning that debugging of postprocessing will disable parallelization
            warnings.warn(
                "Debug output might seem to be redundant due to the parallel processes "
                + "with possibly similar errors in parallel mode. "
                + "Don't be confused!"
            )
        if self.refine.engine == "xtb":
            # Check for f-block elements in forbidden elements
            if self.generate.forbidden_elements:
                if verbosity > 0:
                    f_block_elements = set(range(56, 71)) | set(range(88, 103))
                    lanthanides = set(range(56, 71))
                    if not all(
                        elem in self.generate.forbidden_elements for elem in lanthanides
                    ) or any(
                        elem in f_block_elements
                        for elem in self.generate.element_composition
                    ):
                        warnings.warn(
                            "f-block elements could be within the molecule. "
                            + "xTB does not treat f electrons explicitly. "
                            + "UHF is temporarily set to 0."
                        )

            # Check for super heavy elements in forbidden elements
            super_heavy_elements = set(range(86, 102))
            if self.generate.element_composition and any(
                elem in super_heavy_elements
                for elem in self.generate.element_composition
            ):
                if verbosity > 0:
                    warnings.warn(
                        "xTB does not treat super heavy elements. Atomic numbers are temporarily reduced by 32 to their lighter homologues and then replaced with the correct atom number."
                    )
                # Check if postprocessing is turned off
                if not self.general.postprocess:
                    if verbosity > 0:
                        warnings.warn(
                            "Postprocessing is turned off. The structure will not be relaxed."
                        )

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
