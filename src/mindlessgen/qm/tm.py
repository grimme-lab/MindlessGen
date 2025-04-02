"""
This module handles all Turbomole-related functionality.
"""

from pathlib import Path
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import TURBOMOLEConfig
from .base import QMMethod


class Turbomole(QMMethod):
    """
    This class handles all interaction with the turbomole external dependency.
    """

    def __init__(
        self,
        jobex_path: str | Path,
        ridft_path: str | Path,
        turbomolecfg: TURBOMOLEConfig,
    ) -> None:
        """
        Initialize the turbomole class.
        """
        if isinstance(jobex_path, str):
            self.jobex_path: Path = Path(jobex_path).resolve()
        elif isinstance(jobex_path, Path):
            self.jobex_path = jobex_path
        else:
            raise TypeError("jobex_path should be a string or a Path object.")

        if isinstance(ridft_path, str):
            self.ridft_path: Path = Path(ridft_path).resolve()
        elif isinstance(ridft_path, Path):
            self.ridft_path = ridft_path
        else:
            raise TypeError("ridft_path should be a string or a Path object.")
        self.cfg = turbomolecfg
        self.tmp_dir = self.__class__.get_temporary_directory()

    def optimize(
        self,
        molecule: Molecule,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
    ) -> Molecule:
        """
        Optimize a molecule using Turbomole.
        """
        # Create a unique temporary directory using TemporaryDirectory context manager
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "turbomole_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_coord_to_file(temp_path / "coord")

            # write the input file
            inputname = "control"
            tm_input = self._gen_input(molecule)
            if verbosity > 1:
                print("Turbomole input file:\n##################")
                print(tm_input)
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(tm_input)

            # Setup the turbomole optimization command including the max number of optimization cycles
            arguments = [f"PARNODES={ncores} {self.jobex_path} -c {max_cycles}"]
            if verbosity > 2:
                print(f"Running command: {' '.join(arguments)}")

            tm_log_out, tm_log_err, return_code = self._run_opt(
                temp_path=temp_path, arguments=arguments
            )
            if verbosity > 2:
                print(tm_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"Turbomole failed with return code {return_code}:\n{tm_log_err}"
                )

            # # revert the coord file to xyz file
            coordfile = Path(temp_path / "coord").resolve()
            optimized_molecule = molecule.copy()
            optimized_molecule.read_xyz_from_coord(coordfile)
            return optimized_molecule

    def singlepoint(self, molecule: Molecule, ncores: int, verbosity: int = 1) -> str:
        """
        Perform a single point calculation using Turbomole.
        """
        # Create a unique temporary directory using TemporaryDirectory context manager
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "turbomole_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molfile = "coord"
            molecule.write_coord_to_file(temp_path / molfile)

            # write the input file
            inputname = "control"
            tm_input = self._gen_input(molecule)
            if verbosity > 1:
                print("Turbomole input file:\n##################")
                print(self._gen_input(molecule))
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(tm_input)

            # set up the turbomole single point calculation command
            run_tm = [f"PARNODES={ncores} {self.ridft_path}"]
            if verbosity > 2:
                print(f"Running command: {' '.join(run_tm)}")

            tm_log_out, tm_log_err, return_code = self._run(
                temp_path=temp_path, arguments=run_tm
            )
            if verbosity > 2:
                print(tm_log_out)
            if return_code != 0:
                raise RuntimeError(
                    f"Turbomole failed with return code {return_code}:\n{tm_log_err}"
                )

            return tm_log_out

    def check_gap(
        self, molecule: Molecule, ncores: int, threshold: float, verbosity: int = 1
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
            tm_out = sp.run(
                arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
                shell=True,  # shell=True is necessary for inserting the `PARNODES=xx` command, unfortunately.
            )
            # get the output of the turbomole calculation (of both stdout and stderr)
            tm_log_out = tm_out.stdout.decode("utf8", errors="replace")
            tm_log_err = tm_out.stderr.decode("utf8", errors="replace")

            if (
                "ridft : all done" not in tm_log_out
                or "ridft did not converge!" in tm_log_out
            ):
                raise sp.CalledProcessError(
                    1,
                    str(self.ridft_path),
                    tm_log_out.encode("utf8"),
                    tm_log_err.encode("utf8"),
                )
            return tm_log_out, tm_log_err, 0
        except sp.CalledProcessError as e:
            tm_log_out = e.stdout.decode("utf8", errors="replace")
            tm_log_err = e.stderr.decode("utf8", errors="replace")
            return tm_log_out, tm_log_err, e.returncode

    def _run_opt(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run turbomole optimization with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to turbomole.

        Returns:
        tuple[str, str, int]: The output of the turbomole calculation (stdout and stderr)
                              and the return code
        """
        output_file = temp_path / "job.last"
        try:
            tm_out = sp.run(
                arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
                shell=True,  # shell=True is necessary for inserting the `PARNODES=xx` command, unfortunately.
            )
            # Read the job-last file to get the output of the calculation
            if output_file.exists():
                with open(output_file, encoding="utf-8") as file:
                    tm_log_out = file.read()
                tm_log_err = tm_out.stderr.decode("utf8", errors="replace")
            else:
                raise FileNotFoundError(f"Output file {output_file} not found.")

            if "ridft : all done" not in tm_log_out:
                raise sp.CalledProcessError(
                    1,
                    str(self.jobex_path),
                    tm_log_out,
                    tm_log_err.encode("utf8"),
                )
            return tm_log_out, tm_log_err, 0
        except sp.CalledProcessError as e:
            if output_file.exists():
                with open(output_file, encoding="utf-8", errors="replace") as file:
                    tm_log_out = file.read()
            else:
                tm_log_out = e.stdout.decode(
                    "utf8", errors="replace"
                )  # To prevent the program from crashing a different output is used.
            tm_log_err = e.stderr.decode("utf8", errors="replace")
            return tm_log_out, tm_log_err, e.returncode

    def _gen_input(
        self,
        molecule: Molecule,
    ) -> str:
        """
        Generate a default input file for Turbomole.
        """
        tm_input = "$coord file=coord\n"
        tm_input += f"$charge={molecule.charge} unpaired={molecule.uhf}\n"
        tm_input += "$atoms\n"
        tm_input += f"   basis={self.cfg.basis}\n"
        tm_input += "$dft\n"
        tm_input += f"   functional {self.cfg.functional}\n"
        tm_input += "   gridsize m4\n"
        tm_input += "$rij\n"
        if self.cfg.functional.endswith("-v"):
            tm_input += "$doscnl\n"
        else:
            tm_input += "$disp4\n"
        tm_input += f"$scfiterlimit {self.cfg.scf_cycles}\n"
        tm_input += "$scfconv 7\n"
        tm_input += "$energy file=energy\n"
        tm_input += "$grad file=gradient\n"
        tm_input += "$end"
        return tm_input


# TODO: 1. Convert this to a @staticmethod of Class turbomole
#       2. Rename to `get_method` or similar to enable an abstract interface
#       3. Add the renamed method to the ABC `QMMethod`
#       4. In `main.py`: Remove the passing of the path finder functions as arguments
#          and remove the boiler plate code to make it more general.
def get_ridft_path(binary_name: str | Path | None = None) -> Path:
    """
    Retrieve the path to the Turbomole 'ridft' binary.

    Returns:
    Path: The absolute path to the 'ridft' binary.
    """
    default_ridft_names: list[str | Path] = ["ridft"]
    # put binary name at the beginning of the list to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_ridft_names
    else:
        binary_names = default_ridft_names
    # Get turbomole path from 'which ridft' command
    for binpath in binary_names:
        which_ridft = shutil.which(binpath)
        if which_ridft:
            ridft_path = Path(which_ridft)
            return ridft_path
    raise ImportError("'ridft' (TURBOMOLE DFT) binary could not be found.")


def get_jobex_path(binary_name: str | Path | None = None) -> Path:
    """
    Retrieve the path to the Turbomole 'jobex' binary.

    Returns:
    Path: The absolute path to the 'jobex' binary.
    """
    default_jobex_names: list[str | Path] = ["jobex"]
    # put binary name at the beginning of the list to prioritize it
    if binary_name is not None:
        binary_names = [binary_name] + default_jobex_names
    else:
        binary_names = default_jobex_names
    # Get turbomole path from 'which jobex' command
    for binpath in binary_names:
        which_jobex = shutil.which(binpath)
        if which_jobex:
            jobex_path = Path(which_jobex).resolve()
            return jobex_path
    raise ImportError("'jobex' (TURBOMOLE) binary could not be found.")
