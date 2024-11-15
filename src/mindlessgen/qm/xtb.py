"""
This module handles all xtb-related functionality.
"""

import subprocess as sp
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import numpy as np
from ..molecules import Molecule
from .base import QMMethod
from ..prog import XTBConfig
from ..molecules.miscellaneous import (
    get_lanthanides,
    get_actinides,
)


class XTB(QMMethod):
    """
    This class handles all xtb-related functionality.
    """

    def __init__(self, path: str | Path, xtb_config: XTBConfig) -> None:
        """
        Initialize the XTB class.
        """
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("xtb_path should be a string or a Path object.")
        self.cfg = xtb_config

    def optimize(
        self, molecule: Molecule, max_cycles: int | None = None, verbosity: int = 1
    ) -> Molecule:
        """
        Optimize a molecule using xtb.
        """
        super_heavy_elements = False
        if np.any(molecule.ati > 85):
            ati_original = molecule.ati.copy()
            super_heavy_elements = True
            molecule.ati[molecule.ati > 85] -= 32
        if np.any(np.isin(molecule.ati, get_lanthanides())) or np.any(
            np.isin(molecule.ati, get_actinides())
        ):
            # Store the original UHF value and set uhf to 0
            # Justification: xTB does not treat f electrons explicitly.
            # The remaining openshell system has to be removed.
            uhf_original = molecule.uhf
            molecule.uhf = 0
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
                f"{self.cfg.level}",
            ]
            if molecule.charge != 0:
                arguments += ["--chrg", str(molecule.charge)]
            if molecule.uhf != 0:
                arguments += ["--uhf", str(molecule.uhf)]
            if max_cycles is not None:
                arguments += ["--cycles", str(max_cycles)]

            if verbosity > 2:
                print(f"Running command: {' '.join(arguments)}")

            xtb_log_out, xtb_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(xtb_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"xtb failed with return code {return_code}:\n{xtb_log_err}"
                )

            # read the optimized molecule
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_file(temp_path / "xtbopt.xyz")
            if np.any(np.isin(molecule.ati, get_lanthanides())) or np.any(
                np.isin(molecule.ati, get_actinides())
            ):
                # Reset the UHF value to the original value before returning the optimized molecule.
                optimized_molecule.uhf = uhf_original
            if super_heavy_elements:
                # Reset the atomic numbers to the original values before returning the optimized molecule.
                optimized_molecule.ati = ati_original  # pylint: disable=E0606
                optimized_molecule.atlist = molecule.atlist
            return optimized_molecule

    def singlepoint(self, molecule: Molecule, verbosity: int = 1) -> str:
        """
        Perform a single-point calculation using xtb.
        """
        super_heavy_elements = False
        if np.any(molecule.ati > 85):
            ati_original = molecule.ati.copy()
            super_heavy_elements = True
            molecule.ati[molecule.ati > 85] -= 32
        if np.any(np.isin(molecule.ati, get_lanthanides())) or np.any(
            np.isin(molecule.ati, get_actinides())
        ):
            # Store the original UHF value and set uhf to 0
            # Justification: xTB does not treat f electrons explicitly.
            # The remaining openshell system has to be removed.
            uhf_original = molecule.uhf
            molecule.uhf = 0

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="xtb_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run xtb
            arguments = [
                "molecule.xyz",
                "--gfn",
                f"{self.cfg.level}",
            ]
            if molecule.charge != 0:
                arguments += ["--chrg", str(molecule.charge)]
            if molecule.uhf != 0:
                arguments += ["--uhf", str(molecule.uhf)]

            if verbosity > 2:
                print(f"Running command: xtb {' '.join(arguments)}")

            xtb_log_out, xtb_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(xtb_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"xtb failed with return code {return_code}:\n{xtb_log_err}"
                )

            if np.any(np.isin(molecule.ati, get_lanthanides())) or np.any(
                np.isin(molecule.ati, get_actinides())
            ):
                molecule.uhf = uhf_original
            if super_heavy_elements:
                # Reset the atomic numbers to the original values before returning the optimized molecule.
                molecule.ati = ati_original  # pylint: disable=E0606
            return xtb_log_out

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
        if verbosity > 1:
            print(f"xTB HOMO-LUMO gap: {hlgap:5f}")

        if hlgap > threshold:
            return True
        else:
            return False

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
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
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the xtb calculation (of both stdout and stderr)
            xtb_log_out = xtb_out.stdout.decode("utf8", errors="replace")
            xtb_log_err = xtb_out.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, 0
        except sp.CalledProcessError as e:
            xtb_log_out = e.stdout.decode("utf8", errors="replace")
            xtb_log_err = e.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, e.returncode


# TODO: 1. Convert this to a @staticmethod of Class XTB
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
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
