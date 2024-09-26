"""
This module contains all interactions with the GP3-xTB binary
for next-gen tight-binding calculations.
"""

import subprocess as sp
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from .base import QMMethod


class GP3(QMMethod):
    """
    This class handles all interaction with the GP3 external dependency.
    """

    def __init__(self, path: str | Path) -> None:
        """
        Initialize the GP3 class.
        """
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("gp3_path should be a string or a Path object.")

    def singlepoint(self, molecule: Molecule, verbosity: int = 1) -> str:
        """
        Perform a single-point calculation using GP3-xTB.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="gp3_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run gp3
            arguments = [
                "-c",
                "molecule.xyz",
            ]
            # dump molecule.charge and molecule.uhf to a file
            if molecule.charge != 0:
                with open(temp_path / ".CHRG", "w", encoding="utf-8") as f:
                    f.write(str(molecule.charge))
            if molecule.uhf != 0:
                with open(temp_path / ".UHF", "w", encoding="utf-8") as f:
                    f.write(str(molecule.uhf))

            if verbosity > 2:
                print(f"Running command: gp3 {' '.join(arguments)}")

            gp3_log_out, gp3_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(gp3_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"GP3-xTB failed with return code {return_code}:\n{gp3_log_err}"
                )

            return gp3_log_out

    def check_gap(
        self, molecule: Molecule, threshold: float = 0.5, verbosity: int = 1
    ) -> bool:
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
            gp3_out = self.singlepoint(molecule)
        except RuntimeError as e:
            raise RuntimeError("Single point calculation failed.") from e

        # Parse the output to get the gap
        hlgap = None
        for line in gp3_out.split("\n"):
            if "gap (eV)" in line and "dE" not in line:
                # check if "alpha->alpha" is present in the same line
                # then, the line looks as follows:
                # gap (eV)  alpha->alpha  :        7.02263204
                if "alpha->alpha" in line:
                    hlgap = float(line.split()[4])
                    break
                # otherwise, the line looks as follows:
                # gap (eV)                :        9.06694119
                hlgap = float(line.split()[3])
                break
        if hlgap is None:
            raise ValueError("GP3-xTB gap not determined.")

        if verbosity > 1:
            print(f"GP3-xTB HOMO-LUMO gap: {hlgap:5f}")

        return hlgap > threshold

    def optimize(
        self, molecule: Molecule, max_cycles: int | None = None, verbosity: int = 1
    ) -> Molecule:
        """
        Optimize a molecule using GP3-xTB.
        """
        raise NotImplementedError("Optimization is not yet implemented for GP3-xTB.")

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run GP3-xTB with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to GP3-xTB.

        Returns:
        tuple[str, str, int]: The output of the GP3-xTB calculation (stdout and stderr)
                              and the return code
        """
        try:
            gp3_out = sp.run(
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the GP3-xTB calculation (of both stdout and stderr)
            gp3_log_out = gp3_out.stdout.decode("utf8")
            gp3_log_err = gp3_out.stderr.decode("utf8")
            if (
                "no SCF convergence" in gp3_log_out
                or "nuclear repulsion" not in gp3_log_out
            ):
                raise sp.CalledProcessError(
                    1,
                    str(self.path),
                    gp3_log_out.encode("utf8"),
                    gp3_log_err.encode("utf8"),
                )
            return gp3_log_out, gp3_log_err, 0
        except sp.CalledProcessError as e:
            gp3_log_out = e.stdout.decode("utf8")
            gp3_log_err = e.stderr.decode("utf8")
            return gp3_log_out, gp3_log_err, e.returncode


# TODO: 1. Convert this to a @staticmethod of Class GP3
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
def get_gp3_path(binary_name: str | Path | None = None) -> Path:
    """
    Get the path to the GP3 binary based on different possible names
    that are searched for in the PATH.
    """
    default_gp3_names: list[str | Path] = ["gp3", "gp3_dev"]
    # put binary name at the beginning of the lixt to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_gp3_names
    else:
        binary_names = default_gp3_names
    # Get gp3 path from 'which gp3' command
    for binpath in binary_names:
        which_gp3 = shutil.which(binpath)
        if which_gp3:
            gp3_path = Path(which_gp3).resolve()
            return gp3_path
    raise ImportError("'gp3' binary could not be found.")
