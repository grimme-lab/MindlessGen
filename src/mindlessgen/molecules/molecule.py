"""
Molecule class.
"""

from __future__ import annotations
import copy
from pathlib import Path
import hashlib
import numpy as np
from ..__version__ import __version__


PSE: dict[int, str] = {
    0: "X",
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
PSE_NUMBERS: dict[str, int] = {k.lower(): v for v, k in PSE.items()}
PSE_SYMBOLS: dict[int, str] = {v: k.lower() for v, k in PSE.items()}


class Molecule:
    """
    A class representing a molecule.
    """

    def __init__(self, name: str = ""):
        """
        Initialize a molecule with a name and an optional number of atoms.

        :param name: The name of the molecule.
        :param num_atoms: The initial number of atoms in the molecule (default is 0).
        """
        self._name = name

        self._num_atoms: int | None = None
        self._charge: int | None = None
        self._uhf: int | None = None
        self._atlist: np.ndarray = np.array([], dtype=int)
        self._xyz: np.ndarray = np.array([], dtype=float)
        self._ati: np.ndarray = np.array([], dtype=int)

        self.rng = np.random.default_rng()

    def __str__(self) -> str:
        """
        Return a user-friendly string representation of the molecule.
        """
        returnstr: str = ""
        first_line = True
        if self._name:
            if not first_line:
                returnstr += "\n"
            returnstr += f"Molecule: {self.name}"
            first_line = False
        if self._num_atoms is not None:
            if not first_line:
                returnstr += "\n"
            returnstr += f"# atoms: {self.num_atoms}"
            first_line = False
        if self._charge is not None:
            if not first_line:
                returnstr += "\n"
            returnstr += f"total charge: {self.charge}"
            first_line = False
        if self._uhf is not None:
            if not first_line:
                returnstr += "\n"
            returnstr += f"# unpaired electrons: {self.uhf}"
            first_line = False
        if self._atlist.size:
            if not first_line:
                returnstr += "\n"
            returnstr += f"atomic numbers: {self.atlist}\n"
            returnstr += f"sum formula: {self.sum_formula()}"
            first_line = False
        if self._xyz.size:
            if not first_line:
                returnstr += "\n"
            returnstr += f"atomic coordinates:\n{self.xyz}"
            first_line = False
        if self._ati.size:
            if not first_line:
                returnstr += "\n"
            returnstr += f"atomic number per index: {self._ati}"
            first_line = False
        return returnstr

    def __repr__(self) -> str:
        """
        Return an unambiguous string representation of the molecule.
        """
        # print in this order: name, num_atoms, charge, uhf, xyz
        returnstr = (
            f"Molecule(name={self._name}, "
            + f"num_atoms={self._num_atoms}, "
            + f"charge={self._charge}, "
            + f"uhf={self._uhf}, "
            + f"ati={self._ati}, "
            + f"atlist={self._atlist}, "
            + f"xyz={self._xyz}), "
            + f"sum_formula: {self.sum_formula()}"
        )
        return returnstr

    @staticmethod
    def read_mol_from_file(file: str | Path) -> Molecule:
        """
        Read the XYZ coordinates and the charge of the molecule from a file.
        Thereby, generate a completely new molecule object from scratch.

        Can be called like this:
            from molecule import Molecule
            # Call the static method using the class name
            xyz_file = "example_molecule.xyz"
            molecule_instance = Molecule.read_mol_from_file(xyz_file)
            # Now you can use the molecule_instance as needed
            print(molecule_instance.name)

        The layout of the file is as follows:
        ```
        num_atoms
        'Generated by mindlessgen-v{__version__}'
        <symbol 1> <x1> <y1> <z1>
        <symbol 2> <x2> <y2> <z2>
        ...
        ```

        :param file: The XYZ file to read from.
        :return: A new instance of Molecule with the read data.
        """
        molecule = Molecule()
        if isinstance(file, str):
            file_path = Path(file).resolve()
        elif isinstance(file, Path):
            file_path = file.resolve()
        else:
            raise TypeError("String or Path expected.")
        molecule.read_xyz_from_file(file_path)
        if file_path.with_suffix(".CHRG").exists():
            molecule.read_charge_from_file(file_path.with_suffix(".CHRG"))
        else:
            molecule.charge = 0
        if file_path.with_suffix(".UHF").exists():
            molecule.read_uhf_from_file(file_path.with_suffix(".UHF"))
        else:
            molecule.uhf = 0
        molecule.name = file_path.stem
        return molecule

    @property
    def name(self) -> str:
        """
        Get the name of the molecule.

        :return: The name of the molecule.
        """
        return self._name

    @name.setter
    def name(self, value: str):
        """
        Set the name of the molecule.

        :param value: The name to set.
        :raise TypeError: If the value is not a string.
        """
        if not isinstance(value, str):
            raise TypeError("String expected.")

        self._name = value

    @property
    def num_atoms(self) -> int:
        """
        Get the number of atoms in the molecule.

        :return: The number of atoms in the molecule.
        """
        if self._num_atoms is not None:
            return self._num_atoms
        else:
            if self._atlist.size:
                self.num_atoms = np.sum(self._atlist)
            elif self._xyz.size:
                self.num_atoms = len(self._xyz)
            elif self._ati.size:
                self.num_atoms = len(self._ati)
            if self._num_atoms:
                return self._num_atoms
            else:
                raise ValueError("Number of atoms not present and could not be set.")

    @num_atoms.setter
    # can be either int or numpy.int64
    def num_atoms(self, value: int | np.int64):
        """
        Set the number of atoms in the molecule.

        :param value: The number of atoms to set.
        :raise TypeError: If the value is not an integer.
        :raise ValueError: If the number of atoms is negative.
        """
        if not isinstance(value, int) and not isinstance(value, np.int64):
            raise TypeError("Integer expected.")
        if value < 0:
            raise ValueError("Number of atoms cannot be negative.")

        self._num_atoms = int(value)

    @property
    def charge(self) -> int:
        """
        Get the charge of the molecule.

        :return: The charge of the molecule.
        """
        if self._charge is not None:
            return self._charge
        else:
            raise ValueError("Charge not set.")

    @charge.setter
    def charge(self, value: int | float):
        """
        Set the charge of the molecule.

        :param value: The charge to set.
        :raise TypeError: If the value is not an integer.
        """
        try:
            value = int(value)
        except ValueError as e:
            raise TypeError("Integer expected.") from e

        self._charge = value

    @property
    def uhf(self) -> int:
        """
        Get the UHF of the molecule.

        :return: The UHF of the molecule.
        """
        if self._uhf is not None:
            return self._uhf
        else:
            raise ValueError("UHF not set.")

    @uhf.setter
    def uhf(self, value: int | float):
        """
        Set the UHF of the molecule.

        :param value: The UHF to set.
        :raise TypeError: If the value is not an integer.
        """
        try:
            value = int(value)
        except ValueError as e:
            raise TypeError("Integer expected.") from e

        if value < 0:
            raise ValueError("Number of unpaired electrons cannot be negative.")

        self._uhf = value

    @property
    def xyz(self) -> np.ndarray:
        """
        Get the XYZ coordinates of the molecule.

        :return: The XYZ coordinates of the molecule.
        """
        return self._xyz

    @xyz.setter
    def xyz(self, value: np.ndarray):
        """
        Set the XYZ coordinates of the molecule.

        :param value: The XYZ coordinates to set.
        :raise TypeError: If the value is not a numpy array.
        :raise ValueError: If the shape of the array is not (num_atoms, 3).
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Numpy array expected.")
        if value.shape[1] != 3:
            raise ValueError("Shape of array must be (num_atoms, 3).")

        self._xyz = value

    @property
    def ati(self) -> np.ndarray:
        """
        Get the atomic number per index of the molecule.

        :return: The atomic number per index of the molecule.
        """
        return self._ati

    @ati.setter
    def ati(self, value: np.ndarray):
        """
        Set the atomic number per atom in the molecule.

        :param value: The atomic number per index to set.
        :raise TypeError: If the value is not a numpy array.
        :raise ValueError: If the shape of the array is not (num_atoms,).
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Numpy array expected.")
        # check if array has the right shape
        if value.ndim != 1:
            raise ValueError("Array must have one dimension.")
        # if num_atoms is set, check if the array has the right length
        if self._num_atoms:
            if value.shape[0] != self._num_atoms:
                raise ValueError("Shape of array must be (num_atoms,).")

        self._ati = value

    @property
    def atlist(self) -> np.ndarray:
        """
        Get the initial array with the len 103 (number of elements in the periodic table)
        with the number of atoms of each element.

        :return: The array with the number of atoms of each element.
        """
        return self._atlist

    @atlist.setter
    def atlist(self, value: np.ndarray):
        """
        Set the array with the number of atoms of each element.

        :param value: The atomic numbers to set.
        :raise TypeError: If the value is not a numpy array.
        :raise ValueError: If the shape of the array is not (num_atoms,).
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Numpy array expected.")
        # check if array has the right shape
        if value.ndim != 1:
            raise ValueError("Array must have one dimension.")

        self._atlist = value

    def get_xyz_str(self) -> str:
        """
        Obtain a string with the full XYZ file information of the molecule.
        """
        xyz_str = f"{self.num_atoms}\n"
        try:
            commentline = f"Total charge: {self.charge} ; "
        except ValueError:
            commentline = ""
        try:
            commentline = commentline + f"Unpaired electrons: {self.uhf} ; "
        except ValueError:
            pass
        commentline = commentline + f"Generated by mindlessgen-v{__version__}\n"
        xyz_str += commentline
        for i in range(self.num_atoms):
            xyz_str += (
                f"{PSE[self.ati[i]+1]:<5} "
                + f"{self.xyz[i, 0]:>12.7f} "
                + f"{self.xyz[i, 1]:>12.7f} "
                + f"{self.xyz[i, 2]:>12.7f}\n"
            )
        return xyz_str

    def write_xyz_to_file(self, filename: str | Path | None = None):
        """
        Write the XYZ coordinates of the molecule to a file.

        The layout of the file is as follows:
        ```
        num_atoms
        'Generated by mindlessgen-v{__version__}'
        <symbol 1> <x1> <y1> <z1>
        <symbol 2> <x2> <y2> <z2>
        ...
        ```

        :param filename: The name of the file to write to.
        """
        # raise an error if the number of atoms is not set
        if self._num_atoms is None:
            raise ValueError("Number of atoms not set.")
        if not self._ati.size:
            raise ValueError("Atomic numbers not set.")
        if not self._xyz.size:
            raise ValueError("Atomic coordinates not set.")

        if filename:
            if not isinstance(filename, Path):
                filename = Path(filename).resolve()
        else:
            filename = Path("mlm_" + self.name + ".xyz").resolve()

        with open(filename, "w", encoding="utf8") as f:
            f.write(self.get_xyz_str())
        # if the charge is set, write it to a '.CHRG' file
        if self._charge is not None and self._charge != 0:
            with open(filename.with_suffix(".CHRG"), "w", encoding="utf8") as f:
                f.write(f"{self.charge}\n")
        # if the UHF is set, write it to a '.UHF' file
        if self._uhf is not None and self._uhf > 0:
            with open(filename.with_suffix(".UHF"), "w", encoding="utf8") as f:
                f.write(f"{self.uhf}\n")

    def read_xyz_from_file(self, filename: str | Path) -> None:
        """
        Read the XYZ coordinates of the molecule from a file.

        The layout of the file is as follows:
        ```
        num_atoms
        'Generated by mindlessgen-v{__version__}'
        <symbol 1> <x1> <y1> <z1>
        <symbol 2> <x2> <y2> <z2>
        ...
        ```

        :param filename: The name of the file to read from.
        """
        with open(filename, encoding="utf8") as f:
            lines = f.readlines()
            # read the number of atoms
            self.num_atoms = int(lines[0])
            # read the atomic coordinates
            self.xyz = np.zeros((self.num_atoms, 3))
            self.ati = np.zeros(self.num_atoms, dtype=int)
            self.atlist = np.zeros(103, dtype=int)
            for i in range(self.num_atoms):
                line = lines[i + 2].split()
                self.ati[i] = PSE_NUMBERS[line[0].lower()] - 1
                self.xyz[i, 0] = float(line[1])
                self.xyz[i, 1] = float(line[2])
                self.xyz[i, 2] = float(line[3])
                self.atlist[self.ati[i]] += 1

    def read_charge_from_file(self, filename: str | Path):
        """
        Read the charge of the molecule from a file.

        The layout of the file is as follows:
        ```
        charge
        ```

        :param filename: The name of the file to read from.
        """
        with open(filename, encoding="utf8") as f:
            self.charge = int(f.readline())

    def read_uhf_from_file(self, filename: str | Path):
        """
        Read the UHF of the molecule from a file.

        The layout of the file is as follows:
        ```
        uhf
        ```

        :param filename: The name of the file to read from.
        """
        with open(filename, encoding="utf8") as f:
            self.uhf = int(f.readline())

    def sum_formula(self) -> str:
        """
        Get the sum formula of the molecule.
        """
        if not self._atlist.size:
            raise ValueError("Atomic numbers not set.")
        sumformula = ""
        # begin with C, H, N, O (i.e., 6, 1, 7, 8)
        for i in [5, 0, 6, 7]:
            if self._atlist[i] > 0:
                sumformula += PSE[i + 1] + str(self.atlist[i])
        # Go through all entries of self._ati that are not zero
        for elem, count in enumerate(self.atlist):
            if elem not in [5, 0, 6, 7] and count > 0:
                sumformula += PSE[elem + 1] + str(self.atlist[elem])
        return sumformula

    def copy(self) -> Molecule:
        """
        Create a deep copy of the molecule instance.

        :return: A new instance of Molecule that is a deep copy of the current instance.
        """
        # Create a new instance of Molecule
        new_molecule = Molecule(self._name)

        # Deep copy all attributes
        if self._num_atoms is not None:
            new_molecule.num_atoms = copy.deepcopy(self.num_atoms)
        if self._charge is not None:
            new_molecule.charge = copy.deepcopy(self.charge)
        if self._uhf is not None:
            new_molecule.uhf = copy.deepcopy(self.uhf)
        if self._atlist.size:
            new_molecule.atlist = copy.deepcopy(self.atlist)
        if self._xyz.size:
            new_molecule.xyz = copy.deepcopy(self.xyz)
        if self._ati.size:
            new_molecule.ati = copy.deepcopy(self.ati)

        return new_molecule

    def set_name_from_formula(self) -> None:
        """
        Get the name of the molecule from its sum formula.

        :Arguments: None

        :Returns: None
        """

        molname = self.sum_formula()
        # add a random hash to the name
        hashname = hashlib.sha256(self.rng.bytes(32)).hexdigest()[:6]
        self.name = f"{molname}_{hashname}"


def ati_to_atlist(ati: np.ndarray) -> np.ndarray:
    """
    Convert the atomic number per index to the array with the number of atoms of each element.

    :param ati: The atomic number per index.
    :return: The array with the number of atoms of each element.
    """
    atlist = np.zeros(103, dtype=int)
    for atomtype in ati:
        atlist[atomtype] += 1
    return atlist


def atlist_to_ati(atlist: np.ndarray) -> np.ndarray:
    """
    Convert the array with the number of atoms of each element to the atomic number per index.

    :param atlist: The array with the number of atoms of each element.
    :return: The atomic number per index.
    """
    ati = np.array([], dtype=int)
    for i, num in enumerate(atlist):
        ati = np.append(ati, np.full(shape=num, fill_value=i))
    return ati
