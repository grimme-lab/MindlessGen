"""
This module contains all molecule-related functionality.
"""

from .molecule import Molecule
from .generate_molecule import (
    generate_random_molecule,
    generate_coordinates,
    generate_atom_list,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
    get_alkali_metals,
    get_alkaline_earth_metals,
    check_distances,
)
from .postprocess import iterative_optimization, detect_fragments
from .miscellaneous import set_random_charge

__all__ = [
    "Molecule",
    "generate_random_molecule",
    "generate_coordinates",
    "generate_atom_list",
    "iterative_optimization",
    "detect_fragments",
    "set_random_charge",
    "check_distances",
    "get_three_d_metals",
    "get_four_d_metals",
    "get_five_d_metals",
    "get_lanthanides",
    "get_alkali_metals",
    "get_alkaline_earth_metals",
]
