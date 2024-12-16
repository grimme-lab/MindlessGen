"""
This module contains the Command-Line Interface.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from ..__version__ import __version__


def cli_parser(argv: Sequence[str] | None = None) -> dict:
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    ### General arguments ###
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Input file.",
        required=False,
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        default=False,
        required=False,
        help="Print the configuration and exit.",
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        required=False,
        choices=[-1, 0, 1, 2, 3],
        help="Verbosity level "
        + "(-1 (silent), "
        + "0 (very basic), "
        + "1 (default), "
        + "2 (verbose), or "
        + "3 (super verbose)).",
    )
    parser.add_argument(
        "-P",
        "--parallel",
        type=int,
        required=False,
        help="Number of parallel processes to run.",
    )
    parser.add_argument(
        "-mc",
        "--max-cycles",
        type=int,
        required=False,
        help="Maximum number of optimization cycles.",
    )
    parser.add_argument(
        "--num-molecules",
        type=int,
        required=False,
        help="Number of molecules to generate.",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        default=None,
        required=False,
        help="Postprocess the output.",
    )
    parser.add_argument(
        "--write-xyz",
        action="store_true",
        default=None,
        required=False,
        help="Write the molecules to xyz files.",
    )
    parser.add_argument(
        "--no-write-xyz",
        action="store_false",
        dest="write_xyz",
        required=False,
        help="Do not write the molecules to xyz files.",
    )
    parser.add_argument(
        "--scale-fragment-detection",
        type=float,
        required=False,
        help="Scaling factor for the fragment detection based on the van der Waals radii.",
    )
    parser.add_argument(
        "--scale-minimal-distance",
        type=float,
        required=False,
        help="Minimum atom distance scaling factor.",
    )

    ### Molecule generation arguments ###
    parser.add_argument(
        "--min-num-atoms",
        type=int,
        required=False,
        help="Minimum number of atoms in a molecule.",
    )
    parser.add_argument(
        "--max-num-atoms",
        type=int,
        required=False,
        help="Maximum number of atoms in a molecule.",
    )
    parser.add_argument(
        "--init-coord-scaling",
        type=float,
        required=False,
        help="Initial coordinate scaling factor.",
    )
    parser.add_argument(
        "--increase-scaling-factor",
        type=float,
        required=False,
        help="Factor with which the coordinate scaling factor is increased "
        + "after a failed attempt.",
    )
    parser.add_argument(
        "--element-composition",
        type=str,
        required=False,
        help="Element composition of the molecule.",
    )
    parser.add_argument(
        "--forbidden-elements",
        type=str,
        required=False,
        help="List of forbidden elements.",
    )
    parser.add_argument(
        "--contract_coords",
        type=bool,
        required=False,
        help="Contract the coordinates of the molecule after the coordinats generation.",
    )
    parser.add_argument(
        "--molecular-charge",
        type=str,
        required=False,
        help="Define the charge of the molecule.",
    )
    parser.add_argument(
        "--fixed-composition",
        type=bool,
        required=False,
        help="Fix the element composition of the molecule.",
    )

    ### Refinement arguments ###
    parser.add_argument(
        "--refine-engine",
        type=str,
        required=False,
        choices=["xtb", "orca"],
        help="QM engine to use for refinement.",
    )
    parser.add_argument(
        "--max-frag-cycles",
        type=int,
        required=False,
        help="Maximum number of fragment optimization cycles.",
    )
    parser.add_argument(
        "--refine-hlgap",
        type=float,
        required=False,
        help="Minimum HOMO-LUMO gap required in the end of the refinement process",
    )
    parser.add_argument(
        "--refine-debug",
        action="store_true",
        default=None,
        required=False,
        help="Print debug information during refinement.",
    )

    ### Postprocessing arguments ###
    parser.add_argument(
        "--postprocess-engine",
        type=str,
        required=False,
        choices=["xtb", "orca"],
        help="QM engine to use for postprocessing.",
    )
    parser.add_argument(
        "--postprocess-optimize",
        action="store_true",
        default=None,
        required=False,
        help="Optimize the molecule during post-processing.",
    )
    parser.add_argument(
        "--postprocess-opt-cycles",
        type=int,
        required=False,
        help="Number of optimization cycles in postprocessing.",
    )

    ### xTB specific arguments ###
    parser.add_argument(
        "--xtb-path",
        type=str,
        required=False,
        help="Path to the xTB binary.",
    )
    parser.add_argument(
        "--xtb-level",
        type=int,
        required=False,
        help="Level of theory to use in xTB.",
    )
    parser.add_argument(
        "--postprocess-debug",
        action="store_true",
        default=None,
        required=False,
        help="Print debug information during postprocessing.",
    )

    ### ORCA specific arguments ###
    parser.add_argument(
        "--orca-path",
        type=str,
        required=False,
        help="Path to the ORCA binary.",
    )
    parser.add_argument(
        "--orca-functional",
        type=str,
        required=False,
        help="Functional to use in ORCA.",
    )
    parser.add_argument(
        "--orca-basis",
        type=str,
        required=False,
        help="Basis set to use in ORCA.",
    )
    parser.add_argument(
        "--orca-gridsize",
        type=int,
        required=False,
        help="Solvent to use in ORCA.",
    )
    parser.add_argument(
        "--orca-scf-cycles",
        type=int,
        required=False,
        help="Maximum number of SCF cycles in ORCA.",
    )
    args = parser.parse_args(argv)
    args_dict = vars(args)

    ### TRANSLATE ARGUMENTS TO DICTIONARY ###
    # General arguments
    rev_args_dict: dict[str, dict] = {}
    rev_args_dict["general"] = {
        "config": args_dict["config"],
        "verbosity": args_dict["verbosity"],
        "parallel": args_dict["parallel"],
        "max_cycles": args_dict["max_cycles"],
        "print_config": args_dict["print_config"],
        "num_molecules": args_dict["num_molecules"],
        "postprocess": args_dict["postprocess"],
        "write_xyz": args_dict["write_xyz"],
    }
    # Refinement arguments
    rev_args_dict["refine"] = {
        "engine": args_dict["refine_engine"],
        "max_frag_cycles": args_dict["max_frag_cycles"],
        "hlgap": args_dict["refine_hlgap"],
        "debug": args_dict["refine_debug"],
    }
    # Molecule generation arguments
    rev_args_dict["generate"] = {
        "min_num_atoms": args_dict["min_num_atoms"],
        "max_num_atoms": args_dict["max_num_atoms"],
        "init_coord_scaling": args_dict["init_coord_scaling"],
        "increase_scaling_factor": args_dict["increase_scaling_factor"],
        "element_composition": args_dict["element_composition"],
        "forbidden_elements": args_dict["forbidden_elements"],
        "scale_fragment_detection": args_dict["scale_fragment_detection"],
        "scale_minimal_distance": args_dict["scale_minimal_distance"],
        "contract_coords": args_dict["contract_coords"],
        "molecular_charge": args_dict["molecular_charge"],
        "fixed_composition": args_dict["fixed_composition"],
    }
    # XTB specific arguments
    rev_args_dict["xtb"] = {"xtb_path": args_dict["xtb_path"]}
    rev_args_dict["xtb"]["level"] = args_dict["xtb_level"]
    # ORCA specific arguments
    rev_args_dict["orca"] = {
        "orca_path": args_dict["orca_path"],
        "functional": args_dict["orca_functional"],
        "basis": args_dict["orca_basis"],
        "gridsize": args_dict["orca_gridsize"],
        "scf_cycles": args_dict["orca_scf_cycles"],
    }
    # Postprocessing arguments
    rev_args_dict["postprocess"] = {
        "engine": args_dict["postprocess_engine"],
        "optimize": args_dict["postprocess_optimize"],
        "opt_cycles": args_dict["postprocess_opt_cycles"],
        "debug": args_dict["postprocess_debug"],
    }

    return rev_args_dict
