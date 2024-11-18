"""
This module contains all molecule-related functionality.
"""

from .molecule import Molecule, PSE_NUMBERS, PSE_SYMBOLS, ati_to_atlist, atlist_to_ati
from .generate_molecule import (
    generate_random_molecule,
    generate_coordinates,
    generate_atom_list,
    check_distances,
)
from .refinement import iterative_optimization, detect_fragments
from .postprocess import postprocess_mol
from .miscellaneous import (
    set_random_charge,
    get_three_d_metals,
    get_four_d_metals,
    get_five_d_metals,
    get_lanthanides,
    get_actinides,
    get_alkali_metals,
    get_alkaline_earth_metals,
)

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
    "get_actinides",
    "get_alkali_metals",
    "get_alkaline_earth_metals",
    "PSE_NUMBERS",
    "PSE_SYMBOLS",
    "ati_to_atlist",
    "atlist_to_ati",
    "postprocess_mol",
]
