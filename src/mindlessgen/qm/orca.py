"""
This module handles all ORCA-related functionality.
"""

from collections import defaultdict
from pathlib import Path
import shutil
import subprocess as sp
from tempfile import TemporaryDirectory

from ..molecules import Molecule
from ..prog import DistanceConstraint, ORCAConfig, XTBConfig
from .base import QMMethod
from .xtb import get_xtb_path


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

            inputname = "orca_opt.inp"
            use_xtb_driver = self._should_use_xtb_driver()
            xtb_input = temp_path / "xtb.inp"
            if use_xtb_driver:
                self._write_xtb_input(molecule, xtb_input, inputname)
            orca_input = self._gen_input(
                molecule,
                xyz_filename,
                temp_path,
                ncores,
                True,
                max_cycles,
                use_xtb_driver=use_xtb_driver,
            )
            if verbosity > 1:
                print("ORCA input file:\n##################")
                print(orca_input)
                print("##################")
            with open(temp_path / inputname, "w", encoding="utf8") as f:
                f.write(orca_input)

            # run orca
            if use_xtb_driver:
                orca_log_out, orca_log_err, return_code = self._run_xtb_driver(
                    temp_path=temp_path,
                    geometry_filename=xyz_filename,
                    xcontrol_name=xtb_input.name,
                    ncores=ncores,
                )
            else:
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
            if use_xtb_driver:
                xyzfile = temp_path / "xtbopt.xyz"
                if not xyzfile.exists():
                    raise RuntimeError(
                        "xTB-driven ORCA optimization did not produce 'xtbopt.xyz'."
                    )
            else:
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
        xtb_executable = self._get_xtb_executable()
        arguments = [
            str(xtb_executable),
            geometry_filename,
            "--opt",
        ]
        opt_level = getattr(self.cfg, "optlevel", None)
        if opt_level not in (None, ""):
            arguments.append(str(opt_level))
        arguments.extend(["--orca", "-I", xcontrol_name])
        try:
            xtb_out = sp.run(
                arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            xtb_log_out = xtb_out.stdout.decode("utf8", errors="replace")
            xtb_log_err = xtb_out.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, 0
        except sp.CalledProcessError as e:
            xtb_log_out = e.stdout.decode("utf8", errors="replace")
            xtb_log_err = e.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, e.returncode

    def _get_xtb_executable(self) -> Path:
        """
        Determine the path to the xTB executable for external ORCA optimizations.
        """
        for attr_name in ("xtb_driver_path", "xtb_path"):
            candidate = getattr(self.cfg, attr_name, None)
            if candidate:
                try:
                    return get_xtb_path(candidate)
                except ImportError as exc:
                    raise RuntimeError(
                        f"xTB executable defined via '{attr_name}' could not be found."
                    ) from exc
        try:
            return get_xtb_path(None)
        except ImportError as exc:
            raise RuntimeError(
                "xTB executable not found. Required for constrained ORCA optimizations."
            ) from exc

    def _should_use_xtb_driver(self) -> bool:
        """
        Determine if the xTB external driver should be used (constraints configured).
        """
        return bool(self.xtb_cfg and self.xtb_cfg.distance_constraints)

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
        constraint_lines = self._prepare_distance_constraint_section(molecule)
        lines: list[str] = []
        if constraint_lines:
            lines.append("$constrain")
            if self.xtb_cfg.distance_constraint_force_constant is not None:
                lines.append(
                    f"  force constant= {self.xtb_cfg.distance_constraint_force_constant}"
                )
            lines.extend(constraint_lines)
            lines.append("$end")
        lines.append("$external")
        lines.append(f"  orca input file= {input_file}")
        lines.append(f"  orca bin= {self.path}")
        lines.append("$end")
        xtb_input.write_text("\n".join(lines) + "\n", encoding="utf8")

    def _prepare_distance_constraint_section(self, molecule: Molecule) -> list[str]:
        """
        Convert configured distance constraints to xcontrol instructions.
        """
        if not self.xtb_cfg or not self.xtb_cfg.distance_constraints:
            return []
        element_map: defaultdict[int, list[int]] = defaultdict(list)
        for idx, atomic_number in enumerate(molecule.ati):
            element_map[int(atomic_number)].append(idx)
        constraint_lines: list[str] = []
        for constraint in self.xtb_cfg.distance_constraints:
            self._ensure_constraint_atoms_present(element_map, constraint)
            pairs = self._generate_constraint_pairs(element_map, constraint)
            if not pairs:
                raise RuntimeError(
                    f"No atom pairs found for distance constraint {constraint}."
                )
            for first, second in pairs:
                constraint_lines.append(
                    f"  distance: {first + 1}, {second + 1}, {constraint.distance:.5f}"
                )
        return constraint_lines

    @staticmethod
    def _generate_constraint_pairs(
        element_map: dict[int, list[int]], constraint: DistanceConstraint
    ) -> list[tuple[int, int]]:
        """
        Generate index pairs for the provided constraint.
        """
        atom_a, atom_b = constraint.atomic_numbers
        atom_a_idx = atom_a - 1
        atom_b_idx = atom_b - 1
        indices_a = element_map.get(atom_a_idx, [])
        indices_b = element_map.get(atom_b_idx, [])

        if atom_a == atom_b:
            if len(indices_a) < 2:
                return []
            first, second = sorted(indices_a[:2])
            return [(first, second)]

        if not indices_a or not indices_b:
            return []

        first, second = indices_a[0], indices_b[0]
        if first == second:
            return []
        if first > second:
            first, second = second, first
        return [(first, second)]

    @staticmethod
    def _ensure_constraint_atoms_present(
        element_map: dict[int, list[int]], constraint: DistanceConstraint
    ) -> None:
        """
        Validate that the molecule contains enough atoms for the constraint.
        """
        for atomic_number, required in constraint.required_counts().items():
            idx = atomic_number - 1
            available = len(element_map.get(idx, []))
            if available < required:
                symbol = constraint.symbol_for(atomic_number)
                raise RuntimeError(
                    f"Distance constraint {constraint} requires at least "
                    f"{required} atom(s) of {symbol}, but only {available} present."
                )

    def _gen_input(
        self,
        molecule: Molecule,
        xyzfile: str,
        temp_path: Path,
        ncores: int,
        optimization: bool = False,
        opt_cycles: int | None = None,
        *,
        use_xtb_driver: bool = False,
    ) -> str:
        """
        Generate a default input file for ORCA.
        """
        orca_input = f"! {self.cfg.functional} {self.cfg.basis}\n"
        orca_input += f"! DEFGRID{self.cfg.gridsize}\n"
        orca_input += "! MiniPrint\n"
        orca_input += "! NoTRAH\n"
        if use_xtb_driver:
            orca_input += "! Engrad\n"
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
