"""
This module defines the abstract base class for all QM methods.
"""

from abc import ABC, abstractmethod
from pathlib import Path

from ..molecules import Molecule


class QMMethod(ABC):
    """
    This abstract base class defines the interface for all QM methods.
    """

    def __init__(self):
        pass

    _temp_dir: None | str | Path = None

    # NOTE: Early tests have shown that the parallelization mechanism doesn't support
    # taking the temporary directory as usual class functions/attributes of the QMMethod base class.
    # That's the reason for this seemingly overly complicated implementation.
    # -> For each child class, we access the parent class variable QMMethod._temp_dir
    #    at the time before forking subprocesses,
    #    when accessing the class variable is possible.
    #    The child class then sets the temporary directory (instance variable!)
    #    to the QMMethod._temp_dir
    @classmethod
    def set_temporary_directory(cls, value: str | Path) -> None:
        """
        Set the temporary directory for the QM methods.
        """
        if cls._temp_dir is not None:
            raise ValueError("Parent class variable is already set.")
        if not isinstance(value, (str, Path)):
            raise TypeError("Temporary directory should be a string or a Path object.")
        if isinstance(value, str):
            value = Path(value).resolve()
        elif isinstance(value, Path):
            value = value.resolve()
        # raise error if value is a file
        if value.is_file():
            raise ValueError("Temporary directory should not be a file.")
        value.mkdir(parents=True, exist_ok=True)
        cls._temp_dir = value

    @classmethod
    def get_temporary_directory(cls) -> None | Path:
        """
        Get the temporary directory for the QM methods.
        """
        return cls._temp_dir  # type: ignore[return-value]
        # NOTE: Since a string is always transformed to a Path object in set_temporary_directory,
        #       this is safe to do. The type checker doesn't recognize this, though.

    @abstractmethod
    def optimize(
        self,
        molecule: Molecule,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
    ) -> Molecule:
        """
        Define the optimization process.

        Arguments:
        molecule (Molecule): Molecule to optimize
        ncores (int): Number of cores to use

        Returns:
        Molecule: Optimized molecule
        """

    @abstractmethod
    def singlepoint(self, molecule: Molecule, ncores: int, verbosity: int = 1) -> str:
        """
        Define the single point calculation process.

        Arguments:
        molecule (Molecule): Molecule to calculate
        ncores (int): Number of cores to use
        """

    @abstractmethod
    def check_gap(
        self, molecule: Molecule, ncores: int, threshold: float, verbosity: int = 1
    ) -> bool:
        """
        Check if the HL gap is larger than a given threshold.

        Arguments:
        molecule (Molecule): Molecule to check
        ncores (int): Number of cores to use
        threshold (float): Threshold for the gap
        """

    @abstractmethod
    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Execute the algorithm.

        Arguments:
        temp_path (Path): Path to the temporary directory
        arguments (list[str]): List of arguments to pass to the algorithm

        Returns:
        tuple[str, str, int]: Output, error, and return code
        """
