"""
Molecule class.
"""

from __future__ import annotations

import numpy as np

from ..__version__ import __version__


PSE = {
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
        self._ati: np.ndarray = np.array([], dtype=int)
        self._xyz: np.ndarray = np.array([], dtype=float)

    def __str__(self) -> str:
        """
        Return a user-friendly string representation of the molecule.
        """
        returnstr: str = ""
        if self._name:
            returnstr += f"Molecule: {self.name}\n"
        if self._num_atoms:
            returnstr += f"# atoms: {self.num_atoms}\n"
        if self._charge:
            returnstr += f"total charge: {self.charge}\n"
        if self._uhf:
            returnstr += f"# unpaired electrons: {self.uhf}\n"
        if self._ati.size:
            returnstr += f"atomic numbers: {self.ati}\n"
        if self._xyz.size:
            returnstr += f"atomic coordinates:\n{self.xyz}\n"
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
            + f"xyz={self._xyz})"
        )
        return returnstr

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
        if self._num_atoms is None:
            raise ValueError("Number of atoms not set.")
        return self._num_atoms

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
        if self._charge:
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
        if self._uhf:
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

    def print_xyz(self):
        """
        Print the XYZ coordinates of the molecule.
        """
        print(self._xyz)

    def write_xyz_to_file(self, filename: str):
        """
        Write the XYZ coordinates of the molecule to a file.

        The layout of the file is as follows:
        ```
        num_atoms
        'Generated by mlmgen-v{__version__}'
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

        with open(filename, "w", encoding="utf8") as f:
            f.write(f"{self.num_atoms}\n")
            f.write(f"Generated by mlmgen-v{__version__}\n")
            for i in range(self.num_atoms):
                f.write(
                    f"{PSE[self.ati[i]+1]:<5} "
                    + f"{self.xyz[i, 0]:>12.7f} "
                    + f"{self.xyz[i, 1]:>12.7f} "
                    + f"{self.xyz[i, 2]:>12.7f}\n"
                )

    @property
    def ati(self) -> np.ndarray:
        """
        Get the atomic numbers of the molecule.

        :return: The atomic numbers of the molecule.
        """
        return self._ati

    @ati.setter
    def ati(self, value: np.ndarray):
        """
        Set the atomic numbers of the molecule.

        :param value: The atomic numbers to set.
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

    def sum_formula(self) -> str:
        """
        Get the sum formula of the molecule.
        """

        sumformula = ""
        # begin with C, H, N, O (i.e., 6, 1, 7, 8)
        for i in [5, 0, 6, 7]:
            if self._ati[i] > 0:
                sumformula += PSE[i + 1] + str(self._ati[i])
        # Go through all entries of self._ati that are not zero
        for elem, count in enumerate(self._ati):
            if elem not in [5, 0, 6, 7] and count > 0:
                sumformula += PSE[elem + 1] + str(self._ati[elem])
        return sumformula
