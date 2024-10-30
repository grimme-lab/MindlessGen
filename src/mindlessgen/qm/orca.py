"""
This module handles all ORCA-related functionality.
"""

from pathlib import Path
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import ORCAConfig
from .base import QMMethod


class ORCA(QMMethod):
    """
    This class handles all interaction with the ORCA external dependency.
    """

    def __init__(self, path: str | Path, orcacfg: ORCAConfig) -> None:
        """
        Initialize the ORCA class.
        """
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("orca_path should be a string or a Path object.")
        self.cfg = orcacfg

    def optimize(
        self, molecule: Molecule, max_cycles: int | None = None, verbosity: int = 1
    ) -> Molecule:
        """
        Optimize a molecule using ORCA.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="orca_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(temp_path / "molecule.xyz")

            inputname = "orca_opt.inp"
            orca_input = self._gen_input(molecule, "molecule.xyz", True, max_cycles)
            if verbosity > 1:
                print("ORCA input file:\n##################")
                print(orca_input)
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(orca_input)

            # run orca
            arguments = [
                inputname,
            ]

            orca_log_out, orca_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(orca_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"ORCA failed with return code {return_code}:\n{orca_log_err}"
                )

            # read the optimized molecule from the output file
            xyzfile = Path(temp_path / inputname).resolve().with_suffix(".xyz")
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_file(xyzfile)
            return optimized_molecule

    def singlepoint(self, molecule: Molecule, verbosity: int = 1) -> str:
        """
        Perform a single point calculation using ORCA.
        """
        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="orca_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molfile = "mol.xyz"
            molecule.write_xyz_to_file(temp_path / molfile)

            # write the input file
            inputname = "orca.inp"
            orca_input = self._gen_input(molecule, molfile)
            if verbosity > 1:
                print("ORCA input file:\n##################")
                print(self._gen_input(molecule, molfile))
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(orca_input)

            # run orca
            arguments = [
                inputname,
            ]
            orca_log_out, orca_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(orca_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"ORCA failed with return code {return_code}:\n{orca_log_err}"
                )

            return orca_log_out

    def check_gap(
        self, molecule: Molecule, threshold: float, verbosity: int = 1
    ) -> bool:
        """
        Check if the HL gap is larger than a given threshold.
        """
        raise NotImplementedError("check_gap not implemented for ORCA.")

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run ORCA with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to orca.

        Returns:
        tuple[str, str, int]: The output of the ORCA calculation (stdout and stderr)
                              and the return code
        """
        try:
            orca_out = sp.run(
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the ORCA calculation (of both stdout and stderr)
            orca_log_out = orca_out.stdout.decode("utf8", errors="replace")
            orca_log_err = orca_out.stderr.decode("utf8", errors="replace")
            # check if the output contains "ORCA TERMINATED NORMALLY"
            if "ORCA TERMINATED NORMALLY" not in orca_log_out:
                raise sp.CalledProcessError(
                    1,
                    str(self.path),
                    orca_log_out.encode("utf8"),
                    orca_log_err.encode("utf8"),
                )
            return orca_log_out, orca_log_err, 0
        except sp.CalledProcessError as e:
            orca_log_out = e.stdout.decode("utf8", errors="replace")
            orca_log_err = e.stderr.decode("utf8", errors="replace")
            return orca_log_out, orca_log_err, e.returncode

    def _gen_input(
        self,
        molecule: Molecule,
        xyzfile: str,
        optimization: bool = False,
        opt_cycles: int | None = None,
    ) -> str:
        """
        Generate a default input file for ORCA.
        """
        orca_input = f"! {self.cfg.functional} {self.cfg.basis}\n"
        orca_input += f"! DEFGRID{self.cfg.gridsize}\n"
        orca_input += "! NoTRAH NoSOSCF SlowConv\n"
        # "! AutoAux" keyword for super-heavy elements as def2/J ends at Rn
        if any(atom >= 86 for atom in molecule.ati):
            orca_input += "! AutoAux\n"
        if optimization:
            orca_input += "! OPT\n"
            if opt_cycles is not None:
                orca_input += f"%geom MaxIter {opt_cycles} end\n"
        orca_input += (
            f"%scf\n\tMaxIter {self.cfg.scf_cycles}\n\tConvergence Medium\nend\n"
        )
        orca_input += "%pal nprocs 1 end\n\n"
        orca_input += f"* xyzfile {molecule.charge} {molecule.uhf + 1} {xyzfile}\n"
        return orca_input


# TODO: 1. Convert this to a @staticmethod of Class ORCA
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
def get_orca_path(binary_name: str | Path | None = None) -> Path:
    """
    Get the path to the orca binary based on different possible names
    that are searched for in the PATH.
    """
    default_orca_names: list[str | Path] = ["orca", "orca_dev"]
    # put binary name at the beginning of the lixt to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_orca_names
    else:
        binary_names = default_orca_names
    # Get ORCA path from 'which orca' command
    for binpath in binary_names:
        which_orca = shutil.which(binpath)
        if which_orca:
            orca_path = Path(which_orca).resolve()
            return orca_path
    raise ImportError("'orca' binary could not be found.")
