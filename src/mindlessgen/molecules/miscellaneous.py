"""
Molecule-related helper tools.
"""

import numpy as np


def set_random_charge(
    ati: np.ndarray,
    verbosity: int = 1,
) -> tuple[int, int]:
    """
    Set the charge of a molecule so that unpaired electrons are avoided.
    """
    # go through all ati and add its own value + 1 to get the number of protons
    nel = 0
    for atom in ati:
        nel += atom + 1
    if verbosity > 1:
        print(f"Number of protons in molecule: {nel}")

    if np.any(np.isin(ati, get_lanthanides() + get_actinides())):
        ### Special mode for lanthanides and actinides
        # -> always high spin
        # -> Divide the molecule into Ln3+/Ac3+ ions and negative "ligands"
        # -> The ligands are the remaining protons are assumed to be low spin

        # TODO: use the extracted functions from misscellaneous.py
        uhf = 0
        charge = 0
        ln_protons = 0
        ac_protons = 0
        for atom in ati:
            if atom in get_lanthanides():
                if atom < 64:
                    uhf += atom - 56
                else:
                    uhf += 70 - atom
                ln_protons += (
                    atom - 3 + 1
                )  # subtract 3 to get the number of protons in the Ln3+ ion
            elif atom in get_actinides():
                if atom < 96:
                    uhf += atom - 88
                else:
                    uhf += 102 - atom
                ac_protons += (
                    atom - 3 + 1
                )  # subtract 3 to get the number of protons in the Ln3+ ion
        ligand_protons = nel - ln_protons - ac_protons
        if verbosity > 2:
            if np.any(np.isin(ati, get_lanthanides())):
                print(f"Number of protons from Ln^3+ ions: {ln_protons}")
            if np.any(np.isin(ati, get_actinides())):
                print(f"Number of protons from Ac^3+ ions: {ac_protons}")
            print(
                f"Number of protons from ligands (assuming negative charge): {ligand_protons}"
            )

        if ligand_protons % 2 == 0:
            charge = 0
        else:
            charge = 1
        return charge, uhf
    else:
        ### Default mode
        iseven = False
        if nel % 2 == 0:
            iseven = True
        # if the number of electrons is even, the charge is -2, 0, or 2
        # if the number of electrons is odd, the charge is -1, 1
        rng = np.random.default_rng()
        randint = rng.random()
        if iseven:
            if randint < 1 / 3:
                charge = -2
            elif randint < 2 / 3:
                charge = 0
            else:
                charge = 2
        else:
            if randint < 0.5:
                charge = -1
            else:
                charge = 1
        uhf = 0
        return charge, uhf


def calculate_protons(natoms: np.ndarray) -> int:
    """
    Calculate the number of protons in a molecule from the atom list.
    """
    protons = 0
    for i, atom_count in enumerate(natoms):
        protons += atom_count * (i + 1)
    return protons


def calculate_ligand_electrons(natoms: np.ndarray, nel: int) -> int:
    """
    Calculate the number of ligand electrons in a molecule if lanthanides or actinides are within the molecule.
    """
    f_electrons = sum(
        occurrence
        * (
            ati - 2
        )  # subtract 3 to get the number of electrons for an Ln3+ (Ac3+) ion and add 1 to account for the 0 based indexing.
        for ati, occurrence in enumerate(natoms)
        if (ati in get_lanthanides() or ati in get_actinides())
    )
    ligand_electrons = nel - f_electrons
    return ligand_electrons


def calculate_uhf(atlist: np.ndarray) -> int:
    """
    Calculate the number of unpaired electrons in a molecule.
    """
    uhf = 0
    for ati, occurrence in enumerate(atlist):
        if ati in get_lanthanides():
            if ati < 64:
                uhf += (ati - 56) * occurrence
            else:
                uhf += (70 - ati) * occurrence
        elif ati in get_actinides():
            if ati < 96:
                uhf += (ati - 88) * occurrence
            else:
                uhf += (102 - ati) * occurrence
    return uhf


def get_alkali_metals() -> list[int]:
    """
    Get the atomic numbers of alkali metals.
    """
    alkali = [2, 10, 18, 36, 54]
    return alkali


def get_alkaline_earth_metals() -> list[int]:
    """
    Get the atomic numbers of alkaline earth metals.
    """
    alkaline = [3, 11, 19, 37, 55]
    return alkaline


def get_three_d_metals() -> list[int]:
    """
    Get the atomic numbers of three d metals.
    """
    threedmetals = list(range(20, 30))

    return threedmetals


def get_four_d_metals() -> list[int]:
    """
    Get the atomic numbers of four d metals.
    """

    fourdmetals = list(range(38, 48))
    return fourdmetals


def get_five_d_metals() -> list[int]:
    """
    Get the atomic numbers of five d metals.
    """
    fivedmetals = list(range(71, 80))
    return fivedmetals


def get_lanthanides() -> list[int]:
    """
    Get the atomic numbers of lanthanides.
    """
    lanthanides = list(range(56, 71))
    return lanthanides


def get_actinides() -> list[int]:
    """
    Get the atomic numbers of actinides.
    """
    actinides = list(range(88, 103))
    return actinides
