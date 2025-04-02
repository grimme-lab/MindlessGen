"""
This module contains all interactions with the g-xTB binary
for next-gen tight-binding calculations.
"""

import subprocess as sp
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import GXTBConfig
from .base import QMMethod


class GXTB(QMMethod):
    """
    This class handles all interaction with the g-xTB external dependency.
    """

    def __init__(self, path: str | Path, gxtbcfg: GXTBConfig) -> None:
        """
        Initialize the GXTB class.
        """
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("gxtb_path should be a string or a Path object.")
        self.cfg = gxtbcfg
        self.tmp_dir = self.__class__.get_temporary_directory()

    def singlepoint(self, molecule: Molecule, ncores: int, verbosity: int = 1) -> str:
        """
        Perform a single-point calculation using g-xTB.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "gxtb_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run g-xTB
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
                print(f"Running command: gxtb {' '.join(arguments)}")

            gxtb_log_out, gxtb_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(gxtb_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"g-xTB failed with return code {return_code}:\n{gxtb_log_err}"
                )
            # gp3_output looks like this:
            # [...]
            #   13     -155.03101038        0.00000000        0.00000001       16.45392733   8    F
            #           13  scf iterations
            #           eigenvalues
            # [...]
            # Check for the number of scf iterations
            scf_iterations = 0
            for line in gxtb_log_out.split("\n"):
                if "scf iterations" in line:
                    scf_iterations = int(line.strip().split()[0])
                    break
            if scf_iterations == 0:
                raise RuntimeError("SCF iterations not found in GP3 output.")
            if scf_iterations > self.cfg.scf_cycles:
                raise RuntimeError(
                    f"SCF iterations exceeded limit of {self.cfg.scf_cycles}."
                )

            return gxtb_log_out

    def check_gap(
        self,
        molecule: Molecule,
        ncores: int,
        threshold: float = 0.5,
        verbosity: int = 1,
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
            gxtb_out = self.singlepoint(molecule, ncores)
        except RuntimeError as e:
            raise RuntimeError("Single point calculation failed.") from e

        # Parse the output to get the gap
        hlgap = None
        for line in gxtb_out.split("\n"):
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
            raise ValueError("g-xTB gap not determined.")

        if verbosity > 1:
            print(f"g-xTB HOMO-LUMO gap: {hlgap:5f}")

        return hlgap > threshold

    def optimize(
        self,
        molecule: Molecule,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
    ) -> Molecule:
        """
        Optimize a molecule using g-xTB.
        """
        raise NotImplementedError("Optimization is not yet implemented for g-xTB.")

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run g-xTB with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to g-xTB.

        Returns:
        tuple[str, str, int]: The output of the g-xTB calculation (stdout and stderr)
                              and the return code
        """
        try:
            gxtb_out = sp.run(
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the g-xTB calculation (of both stdout and stderr)
            gxtb_log_out = gxtb_out.stdout.decode("utf8")
            gxtb_log_err = gxtb_out.stderr.decode("utf8")
            if (
                "no SCF convergence" in gxtb_log_out
                or "nuclear repulsion" not in gxtb_log_out
            ):
                raise sp.CalledProcessError(
                    1,
                    str(self.path),
                    gxtb_log_out.encode("utf8"),
                    gxtb_log_err.encode("utf8"),
                )
            return gxtb_log_out, gxtb_log_err, 0
        except sp.CalledProcessError as e:
            gxtb_log_out = e.stdout.decode("utf8")
            gxtb_log_err = e.stderr.decode("utf8")
            return gxtb_log_out, gxtb_log_err, e.returncode


# TODO: 1. Convert this to a @staticmethod of Class GXTB
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
def get_gxtb_path(binary_name: str | Path | None = None) -> Path:
    """
    Get the path to the g-xTB binary based on different possible names
    that are searched for in the PATH.
    """
    default_gxtb_names: list[str | Path] = ["gxtb", "gxtb_dev"]
    # put binary name at the beginning of the lixt to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_gxtb_names
    else:
        binary_names = default_gxtb_names
    # Get g-xTB path from 'which gxtb' command
    for binpath in binary_names:
        which_gxtb = shutil.which(binpath)
        if which_gxtb:
            gxtb_path = Path(which_gxtb).resolve()
            return gxtb_path
    raise ImportError("'gxtb' binary could not be found.")
