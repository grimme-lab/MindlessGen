"""
This module handles all xtb-related functionality.
"""

import subprocess as sp
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from ..molecules import Molecule
from .base import QMMethod


class XTB(QMMethod):
    """
    This class handles all xtb-related functionality.
    """

    def __init__(self, path: str | Path = "xtb", verbosity: int = 1):
        """
        Initialize the XTB class.
        """
        if isinstance(path, str):
            self.xtb_path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.xtb_path = path
        else:
            raise TypeError("xtb_path should be a string or a Path object.")
        self.verbosity = verbosity

    def optimize(self, molecule: Molecule, verbosity: int = 1) -> Molecule:
        """
        Optimize a molecule using xtb.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="xtb_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run xtb
            arguments = [
                "molecule.xyz",
                "--opt",
                "--gfn",
                "2",
            ]
            if self.verbosity > 2:
                print(f"Running command: {' '.join(arguments)}")

            xtb_log_out, xtb_log_err, return_code = self.run(
                temp_path=temp_path, arguments=arguments
            )
            if return_code != 0:
                if verbosity > 2:
                    print(xtb_log_out)
                raise RuntimeError(
                    f"xtb failed with return code {return_code}:\n{xtb_log_err}"
                )

            # read the optimized molecule
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_file(temp_path / "xtbopt.xyz")

            return optimized_molecule

    def singlepoint(self, molecule: Molecule) -> str:
        """
        Optimize a molecule using xtb.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="xtb_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run xtb
            arguments = [
                "molecule.xyz",
                "--gfn",
                "2",
            ]
            if self.verbosity > 1:
                print(f"Running command: {' '.join(arguments)}")

            xtb_log_out, xtb_log_err, return_code = self.run(
                temp_path=temp_path, arguments=arguments
            )
            if return_code != 0:
                raise RuntimeError(
                    f"xtb failed with return code {return_code}:\n{xtb_log_err}"
                )

            return xtb_log_out

    def check_gap(self, molecule: Molecule, threshold: float = 0.5) -> bool:
        """
        Check if the HL gap is larger than a given threshold.

        Arguments:
        molecule (Molecule): Molecule to check
        threshold (float): Threshold for the gap

        Returns:
        bool: True if the gap is larger than the threshold, False otherwise
        """

        # Perform a single point calculation
        try:
            xtb_out = self.singlepoint(molecule)
        except RuntimeError as e:
            raise RuntimeError("Single point calculation failed.") from e

        # Parse the output to get the gap
        hlgap = None
        for line in xtb_out.split("\n"):
            if "HOMO-LUMO GAP" in line:
                hlgap = float(line.split()[3])
                break

        if hlgap is None:
            raise ValueError("HOMO-LUMO gap not determined.")
        if self.verbosity > 1:
            print(f"HOMO-LUMO gap: {hlgap:5f}")

        if hlgap > threshold:
            return True
        else:
            return False

    def run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run xtb with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to xtb.

        Returns:
        tuple[str, str, int]: The output of the xtb calculation (stdout and stderr)
                              and the return code
        """
        non_parallel = ["-P", "1"]
        arguments += non_parallel
        try:
            xtb_out = sp.run(
                [str(self.xtb_path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the xtb calculation (of both stdout and stderr)
            xtb_log_out = xtb_out.stdout.decode("utf8")
            xtb_log_err = xtb_out.stderr.decode("utf8")
            return xtb_log_out, xtb_log_err, 0
        except sp.CalledProcessError as e:
            xtb_log_out = e.stdout.decode("utf8")
            xtb_log_err = e.stderr.decode("utf8")
            return xtb_log_out, xtb_log_err, e.returncode


def get_xtb_path(binary_name: str | Path | None = None) -> Path:
    """
    Get the path to the xtb binary based on different possible names
    that are searched for in the PATH.
    """
    default_xtb_names: list[str | Path] = ["xtb", "xtb_dev"]
    # put binary name at the beginning of the lixt to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_xtb_names
    else:
        binary_names = default_xtb_names
    # Get xtb path from 'which xtb' command
    for binpath in binary_names:
        which_xtb = shutil.which(binpath)
        if which_xtb:
            xtb_path = Path(which_xtb).resolve()
            return xtb_path
    raise ImportError("'xtb' binary could not be found.")
