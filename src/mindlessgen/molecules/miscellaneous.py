"""
Molecule-related helper tools.
"""

import numpy as np


def set_random_charge(ati: np.ndarray, verbosity: int = 1) -> int:
    """
    Set the charge of a molecule so that unpaired electrons are avoided.
    """
    # go through all ati and add its own value + 1 to get the number of protons
    nel = 0
    for atom in ati:
        nel += atom + 1
    if verbosity > 1:
        print(f"Number of protons in molecule: {nel}")
    iseven = False
    if nel % 2 == 0:
        iseven = True
    # if the number of electrons is even, the charge is -2, 0, or 2
    # if the number of electrons is odd, the charge is -1, 1
    randint = np.random.rand()
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
    return charge
