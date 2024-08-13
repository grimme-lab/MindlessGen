"""
This module handles all xtb-related functionality.
"""

import subprocess as sp
from pathlib import Path
import shutil

from tempfile import TemporaryDirectory

from ..molecules import Molecule


class XTB:
    """
    This class handles all xtb-related functionality.
    """

    def __init__(self, xtb_path: str | Path = "xtb", verbosity: int = 1):
        if isinstance(xtb_path, str):
            self.xtb_path: Path = Path(xtb_path).resolve()
        elif isinstance(xtb_path, Path):
            self.xtb_path = xtb_path
        else:
            raise TypeError("xtb_path should be a string or a Path object.")
        self.verbosity = verbosity

    def optimize(self, molecule: Molecule) -> Molecule:
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
            if self.verbosity > 1:
                print(f"Running command: {' '.join(arguments)}")

            xtb_log_out, xtb_log_err, return_code = self.run_xtb(
                temp_path=temp_path, arguments=arguments
            )
            if return_code != 0:
                raise RuntimeError(
                    f"xtb failed with return code {return_code}:\n{xtb_log_err}"
                )

            # read the optimized molecule
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_file(temp_path / "xtbopt.xyz")

            return optimized_molecule

    def run_xtb(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run xtb with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to xtb.

        Returns:
        tuple[str, str, int]: The output of the xtb calculation (stdout and stderr)
                              and the return code
        """
        try:
            xtb_out = sp.run(
                [str(self.xtb_path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the xtb calculation (of both stdout and stderr)
            xtb_log_out = xtb_out.stdout.decode("utf-8")
            xtb_log_err = xtb_out.stderr.decode("utf-8")
            return xtb_log_out, xtb_log_err, 0
        except sp.CalledProcessError as e:
            xtb_log_out = e.stdout.decode("utf-8")
            xtb_log_err = e.stderr.decode("utf-8")
            return xtb_log_out, xtb_log_err, e.returncode


def get_xtb_path(binary_names: list[str]) -> Path | None:
    """
    Get the path to the xtb binary based on different possible names
    that are searched for in the PATH.
    """
    # Get xtb path from 'which xtb' command
    for binary_name in binary_names:
        which_xtb = shutil.which(binary_name)
        if which_xtb:
            xtb_path = Path(which_xtb)
            return xtb_path
    print("xtb not found.")
    return None
