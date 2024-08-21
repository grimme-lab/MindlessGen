from abc import ABC, abstractmethod
from pathlib import Path

from ..molecules import Molecule


class QMMethod(ABC):
    """
    This abstract base class defines the interface for all QM methods.
    """

    @abstractmethod
    def __init__(self, path: str | Path, verbosity: int = 1):
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("xtb_path should be a string or a Path object.")
        self.verbosity = verbosity

    @abstractmethod
    def optimize(
        self, molecule: Molecule, max_cycles: int | None = None, verbosity: int = 1
    ) -> Molecule:
        """
        Define the optimization process.

        Arguments:
        molecule (Molecule): Molecule to optimize

        Returns:
        Molecule: Optimized molecule
        """

    @abstractmethod
    def singlepoint(self, molecule: Molecule, verbosity: int = 1) -> str:
        """
        Define the single point calculation process.

        Arguments:
        molecule (Molecule): Molecule to calculate
        """

    @abstractmethod
    def check_gap(
        self, molecule: Molecule, threshold: float, verbosity: int = 1
    ) -> bool:
        """
        Check if the HL gap is larger than a given threshold.

        Arguments:
        molecule (Molecule): Molecule to check
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
