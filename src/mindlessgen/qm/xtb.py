"""
This module handles all xtb-related functionality.
"""

import subprocess as sp
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from collections import defaultdict
import numpy as np
from ..molecules import Molecule
from .base import QMMethod
from ..prog import XTBConfig, DistanceConstraint
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
        self.tmp_dir = self.__class__.get_temporary_directory()

    def optimize(
        self,
        molecule: Molecule,
        ncores: int,
        max_cycles: int | None = None,
        verbosity: int = 1,
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
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "xtb_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run xtb
            arguments = [
                "molecule.xyz",
                "--opt",
                "--gfn",
                f"{self.cfg.level}",
                "-P",
                f"{ncores}",
            ]
            if molecule.charge != 0:
                arguments += ["--chrg", str(molecule.charge)]
            if molecule.uhf != 0:
                arguments += ["--uhf", str(molecule.uhf)]
            if max_cycles is not None:
                arguments += ["--cycles", str(max_cycles)]
            if self.cfg.distance_constraints:
                if verbosity > 1:
                    print("Preparing distance constraint file...")
                if self._prepare_distance_constraint_file(molecule, temp_path):
                    if verbosity > 1:
                        print("Distance constraint file prepared.")
                    arguments += ["--input", "xtb.inp"]

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

    def singlepoint(self, molecule: Molecule, ncores: int, verbosity: int = 1) -> str:
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
        kwargs_temp_dir: dict[str, str | Path] = {"prefix": "xtb_"}
        if self.tmp_dir is not None:
            kwargs_temp_dir["dir"] = self.tmp_dir
        with TemporaryDirectory(**kwargs_temp_dir) as temp_dir:  # type: ignore[call-overload]
            temp_path = Path(temp_dir).resolve()
            # write the molecule to a temporary file
            molecule.write_xyz_to_file(str(temp_path / "molecule.xyz"))

            # run xtb
            arguments = [
                "molecule.xyz",
                "--gfn",
                f"{self.cfg.level}",
                "-P",
                f"{ncores}",
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
        ncores (int): Number of cores to use
        threshold (float): Threshold for the gap

        Returns:
        bool: True if the gap is larger than the threshold, False otherwise
        """

        # Perform a single point calculation
        try:
            xtb_out = self.singlepoint(molecule, ncores)
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

    def _prepare_distance_constraint_file(
        self, molecule: Molecule, temp_path: Path
    ) -> bool:
        """
        Write an xtb.inp file describing user supplied distance constraints.
        """
        element_map: defaultdict[int, list[int]] = defaultdict(list)
        for idx, atomic_number in enumerate(molecule.ati):
            element_map[int(atomic_number)].append(idx)

        constraint_lines: list[str] = []
        for constraint in self.cfg.distance_constraints:
            self._ensure_constraint_atoms_present(element_map, constraint)
            pairs = self._generate_constraint_pairs(element_map, constraint)
            if not pairs:
                raise RuntimeError(
                    f"No atom pairs found for distance constraint {constraint}."
                )
            for first, second in pairs:
                constraint_lines.append(
                    f" distance: {first + 1}, {second + 1}, {constraint.distance}"
                )

        if not constraint_lines:
            return False

        contents = ["$constrain"]
        if self.cfg.distance_constraint_force_constant is not None:
            contents.append(
                f" force constant= {self.cfg.distance_constraint_force_constant}"
            )
        contents.extend(constraint_lines)
        contents.append("$end")
        (temp_path / "xtb.inp").write_text("\n".join(contents) + "\n", encoding="utf8")
        return True

    @staticmethod
    def _generate_constraint_pairs(
        element_map: dict[int, list[int]], constraint: DistanceConstraint
    ) -> list[tuple[int, int]]:
        """
        Generate the index pair for the provided constraint.
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

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
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
