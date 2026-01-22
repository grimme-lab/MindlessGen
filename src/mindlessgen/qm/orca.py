"""
This module handles all ORCA-related functionality.
"""

from pathlib import Path
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import ORCAConfig, XTBConfig
from .base import QMMethod
from .xtb import XTB, get_xtb_path


class ORCA(QMMethod):
    """
    This class handles all interaction with the ORCA external dependency.
    """

    def __init__(
        self, path: str | Path, orcacfg: ORCAConfig, xtb_config: XTBConfig | None = None
    ) -> None:
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
        self.xtb_cfg = xtb_config
        # must be explicitly initialized in current parallelization implementation
        # as accessing parent class variables might not be possible
        self.tmp_dir = self.__class__.get_temporary_directory()

    def optimize(
        self,
        molecule: Molecule,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
    ) -> Molecule:
        """
        Optimize a molecule using ORCA.
        """

        # Create a unique temporary directory using TemporaryDirectory context manager
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "orca_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            # NOTE: "prefix" and "dir" are valid keyword arguments for TemporaryDirectory
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            xyz_filename = "molecule.xyz"
            molecule.write_xyz_to_file(temp_path / xyz_filename)

            if self.cfg.use_xtb_driver:
                optimized_molecule = self.optimize_xtb_driver(
                    temp_path=temp_path,
                    molecule=molecule,
                    xyz_filename=xyz_filename,
                    ncores=ncores,
                    max_cycles=max_cycles,
                    verbosity=verbosity,
                )
                return optimized_molecule
            inputname = "orca_opt.inp"
            orca_input = self._gen_input(
                molecule,
                xyz_filename,
                temp_path,
                ncores,
                True,
                max_cycles,
            )
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

    def singlepoint(self, molecule: Molecule, ncores: int, verbosity: int = 1) -> str:
        """
        Perform a single point calculation using ORCA.
        """
        # Create a unique temporary directory using TemporaryDirectory context manager
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "orca_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            # NOTE: "prefix" and "dir" (also as Path) are valid keyword arguments
            # for TemporaryDirectory
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molfile = "mol.xyz"
            molecule.write_xyz_to_file(temp_path / molfile)

            # write the input file
            inputname = "orca.inp"
            orca_input = self._gen_input(molecule, molfile, temp_path, ncores)
            if verbosity > 1:
                print("ORCA input file:\n##################")
                print(self._gen_input(molecule, molfile, temp_path, ncores))
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
        self, molecule: Molecule, ncores: int, threshold: float, verbosity: int = 1
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
        _temp_path: Path,
        ncores: int,
        optimization: bool = False,
        opt_cycles: int | None = None,
    ) -> str:
        """
        Generate a default input file for ORCA.
        """
        orca_input = f"! {self.cfg.functional} {self.cfg.basis}\n"
        orca_input += f"! DEFGRID{self.cfg.gridsize}\n"
        orca_input += "! MiniPrint\n"
        orca_input += "! NoTRAH\n"
        # "! AutoAux" keyword for super-heavy elements as def2/J ends at Rn
        if any(atom >= 86 for atom in molecule.ati):
            orca_input += "! AutoAux\n"
        if optimization:
            orca_input += "! OPT\n"
            if opt_cycles is not None:
                orca_input += f"%geom MaxIter {opt_cycles} end\n"
        orca_input += f"%scf\n\tMaxIter {self.cfg.scf_cycles}\n"
        if not optimization:
            orca_input += "\tConvergence Medium\n"
        orca_input += "end\n"
        orca_input += f"%pal nprocs {ncores} end\n\n"
        orca_input += f"* xyzfile {molecule.charge} {molecule.uhf + 1} {xyzfile}\n"
        return orca_input

    def optimize_xtb_driver(
        self,
        temp_path: Path,
        molecule: Molecule,
        xyz_filename: str,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
    ) -> Molecule:
        """
        Optimize a molecule using ORCA through the xTB external driver.
        """

        xtb_input = temp_path / "xtb.inp"
        inputname = "orca_opt.inp"
        self._write_xtb_input(molecule, xtb_input, inputname)
        orca_input = self._gen_input_xtb_driver(
            molecule,
            xyz_filename,
            temp_path,
            ncores,
            True,
            max_cycles,
        )
        if verbosity > 1:
            print("ORCA input file:\n##################")
            print(orca_input)
            print("##################")
            print("XTB input file:\n##################")
            print(xtb_input)
            print("##################")
        with open(temp_path / inputname, "w", encoding="utf8") as f:
            f.write(orca_input)
        # run orca with xTB as a driver
        orca_log_out, orca_log_err, return_code = self._run_xtb_driver(
            temp_path=temp_path,
            geometry_filename=xyz_filename,
            xcontrol_name=xtb_input.name,
            ncores=ncores,
        )
        if verbosity > 2:
            print(orca_log_out)
        if return_code != 0:
            raise RuntimeError(
                f"ORCA failed with return code {return_code}:\n{orca_log_err}"
            )

        # read the optimized molecule from the output file
        xyzfile = temp_path / "xtbopt.xyz"
        if not xyzfile.exists():
            raise RuntimeError(
                "xTB-driven ORCA optimization did not produce 'xtbopt.xyz'."
            )
        optimized_molecule = molecule.copy()
        optimized_molecule.read_xyz_from_file(xyzfile)
        return optimized_molecule

    def _run_xtb_driver(
        self,
        temp_path: Path,
        geometry_filename: str,
        xcontrol_name: str,
        ncores: int,
    ) -> tuple[str, str, int]:
        """
        Run the optimization through the xTB external driver when constraints are requested.
        """
        xtb_executable = get_xtb_path()
        if self.xtb_cfg is None:
            raise RuntimeError(
                "xTB driver requested but no xTB configuration provided."
            )
        xtb_runner = XTB(path=xtb_executable, xtb_config=self.xtb_cfg)
        arguments = [
            geometry_filename,
            "--opt",
        ]
        opt_level = getattr(self.cfg, "optlevel", None)
        if opt_level not in (None, ""):
            arguments.append(str(opt_level))
        arguments.extend(["--orca", "-I", xcontrol_name])
        xtb_log_out, xtb_log_err, returncode = xtb_runner._run(
            temp_path=temp_path, arguments=arguments
        )
        return xtb_log_out, xtb_log_err, returncode

    def _write_xtb_input(
        self, molecule: Molecule, xtb_input: Path, input_file: str
    ) -> None:
        """
        Write the xcontrol file containing constraints and ORCA driver info.
        """
        if not self.xtb_cfg:
            raise RuntimeError(
                "xTB configuration missing but constraints were requested."
            )
        xtb_path = get_xtb_path()
        xtb_writer = XTB(xtb_path, self.xtb_cfg)
        generated = xtb_writer._prepare_distance_constraint_file(
            molecule, xtb_input.parent
        )
        if not generated:
            raise RuntimeError(
                "xTB driver requested but no distance constraints were generated."
            )
        with xtb_input.open("a", encoding="utf8") as handle:
            handle.write("$external\n")
            handle.write(f"  orca input file= {input_file}\n")
            handle.write(f"  orca bin= {self.path}\n")
            handle.write("$end\n")

    def _gen_input_xtb_driver(
        self,
        molecule: Molecule,
        xyzfile: str,
        temp_path: Path,
        ncores: int,
        optimization: bool = False,
        opt_cycles: int | None = None,
    ) -> str:
        """
        Generate a default input file for ORCA.
        """
        orca_input = f"! {self.cfg.functional} {self.cfg.basis}\n"
        orca_input += f"! DEFGRID{self.cfg.gridsize}\n"
        orca_input += "! MiniPrint\n"
        orca_input += "! NoTRAH\n"
        orca_input += "! Engrad\n"
        # "! AutoAux" keyword for super-heavy elements as def2/J ends at Rn
        if any(atom >= 86 for atom in molecule.ati):
            orca_input += "! AutoAux\n"
        orca_input += f"%scf\n\tMaxIter {self.cfg.scf_cycles}\n"
        if not optimization:
            orca_input += "\tConvergence Medium\n"
        orca_input += "end\n"
        orca_input += f"%pal nprocs {ncores} end\n\n"
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
