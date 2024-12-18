"""
This module handles all Turbomole-related functionality.
"""

from pathlib import Path
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import TURBOMOLEConfig as turbomoleConfig
from .base import QMMethod


class Turbomole(QMMethod):
    """
    This class handles all interaction with the turbomole external dependency.
    """

    def __init__(self, path: str | Path, turbomolecfg: turbomoleConfig) -> None:
        """
        Initialize the turbomole class.
        """
        if isinstance(path, str):
            self.path: Path = Path(path).resolve()
        elif isinstance(path, Path):
            self.path = path
        else:
            raise TypeError("turbomole_path should be a string or a Path object.")
        self.cfg = turbomolecfg

    def optimize(
        self, molecule: Molecule, max_cycles: int | None = None, verbosity: int = 1
    ) -> Molecule:
        """
        Optimize a molecule using turbomole.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="turbomole_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(temp_path / "molecule.xyz")

            inputname = "turbomole_opt.inp"
            turbomole_input = self._gen_input(
                molecule, "molecule.xyz", True, max_cycles
            )
            if verbosity > 1:
                print("turbomole input file:\n##################")
                print(turbomole_input)
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(turbomole_input)

            # run turbomole
            arguments = [
                inputname,
            ]

            turbomole_log_out, turbomole_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(turbomole_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"turbomole failed with return code {return_code}:\n{turbomole_log_err}"
                )

            # read the optimized molecule from the output file
            xyzfile = Path(temp_path / inputname).resolve().with_suffix(".xyz")
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_file(xyzfile)
            return optimized_molecule

    def singlepoint(self, molecule: Molecule, verbosity: int = 1) -> str:
        """
        Perform a single point calculation using turbomole.
        """
        # Create a unique temporary directory using TemporaryDirectory context manager
        with TemporaryDirectory(prefix="turbomole_") as temp_dir:
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molfile = "mol.xyz"
            molecule.write_xyz_to_file(temp_path / molfile)

            # converte mol.xyz to tm format
            arguments = [
                "x2t mol.xyz > coord",
            ]
            print("Running", arguments)

            # run cefine to write control file
            arguments = [
                "cefine",
                "-bas",
                f"{self.cfg.basis}",
                "-func",
                f"{self.cfg.functional}",
            ]
            if molecule.charge != 0:
                arguments += ["-chrg", str(molecule.charge)]
            if molecule.uhf != 0:
                arguments += ["-uhf", str(molecule.uhf)]
            # if molecule.symmetry is not None:
            #     arguments += ["-sym", str(molecule.symmetry)]
            print("Running cefine", arguments)

            # run turbomole
            arguments = [
                "ridft > rifdt.out",
            ]
            turbomole_log_out, turbomole_log_err, return_code = self._run(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(turbomole_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"turbomole failed with return code {return_code}:\n{turbomole_log_err}"
                )

            return turbomole_log_out

    def check_gap(
        self, molecule: Molecule, threshold: float, verbosity: int = 1
    ) -> bool:
        """
        Check if the HL gap is larger than a given threshold.
        """
        raise NotImplementedError("check_gap not implemented for turbomole.")

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run turbomole with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to turbomole.

        Returns:
        tuple[str, str, int]: The output of the turbomole calculation (stdout and stderr)
                              and the return code
        """
        try:
            turbomole_out = sp.run(
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the turbomole calculation (of both stdout and stderr)
            turbomole_log_out = turbomole_out.stdout.decode("utf8", errors="replace")
            turbomole_log_err = turbomole_out.stderr.decode("utf8", errors="replace")
            # check if the output contains "turbomole TERMINATED NORMALLY"
            if "all done" not in turbomole_log_out:
                raise sp.CalledProcessError(
                    1,
                    str(self.path),
                    turbomole_log_out.encode("utf8"),
                    turbomole_log_err.encode("utf8"),
                )
            return turbomole_log_out, turbomole_log_err, 0
        except sp.CalledProcessError as e:
            turbomole_log_out = e.stdout.decode("utf8", errors="replace")
            turbomole_log_err = e.stderr.decode("utf8", errors="replace")
            return turbomole_log_out, turbomole_log_err, e.returncode

    def _cefine(self, molecule: Molecule) -> list[str]:
        """
        Refine a molecule using turbomole.
        """
        call = [
            "cefine",
            "-bas",
            f"{self.cfg.basis}",
            "-func",
            f"{self.cfg.functional}",
        ]
        if molecule.charge != 0:
            call += ["-chrg", str(molecule.charge)]
        if molecule.uhf != 0:
            call += ["-uhf", str(molecule.uhf)]
        # if molecule.symmetry is not None:
        #     arguments += ["-sym", str(molecule.symmetry)]

        return call

    def _gen_input(
        self,
        molecule: Molecule,
        xyzfile: str,
        optimization: bool = False,
        opt_cycles: int | None = None,
    ) -> str:
        """
        Generate a default input file for turbomole.
        """
        turbomole_input = f"! {self.cfg.functional} {self.cfg.basis}\n"
        # turbomole_input += f"! DEFGRID{self.cfg.gridsize}\n"
        turbomole_input += "! NoTRAH NoSOSCF SlowConv\n"
        # "! AutoAux" keyword for super-heavy elements as def2/J ends at Rn
        if any(atom >= 86 for atom in molecule.ati):
            turbomole_input += "! AutoAux\n"
        if optimization:
            turbomole_input += "! OPT\n"
            if opt_cycles is not None:
                turbomole_input += f"%geom MaxIter {opt_cycles} end\n"
        turbomole_input += (
            f"%scf\n\tMaxIter {self.cfg.scf_cycles}\n\tConvergence Medium\nend\n"
        )
        turbomole_input += "%pal nprocs 1 end\n\n"
        turbomole_input += f"* xyzfile {molecule.charge} {molecule.uhf + 1} {xyzfile}\n"
        return turbomole_input


# TODO: 1. Convert this to a @staticmethod of Class turbomole
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
def get_turbomole_path(binary_name: str | Path | None = None) -> Path:
    """
    Get the path to the turbomole binary based on different possible names
    that are searched for in the PATH.
    """
    default_turbomole_names: list[str | Path] = ["turbomole", "turbomole_dev"]
    # put binary name at the beginning of the lixt to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_turbomole_names
    else:
        binary_names = default_turbomole_names
    # Get turbomole path from 'which turbomole' command
    for binpath in binary_names:
        which_turbomole = shutil.which(binpath)
        if which_turbomole:
            turbomole_path = Path(which_turbomole).resolve()
            return turbomole_path
    raise ImportError("'turbomole' binary could not be found.")
