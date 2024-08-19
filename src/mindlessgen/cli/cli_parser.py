from __future__ import annotations

import argparse
from collections.abc import Sequence
from ..__version__ import __version__


def cli_parser(argv: Sequence[str] | None = None) -> dict:
    """
    Parse command line arguments.
    """
    # get command line argument
    parser = argparse.ArgumentParser()
    # General arguments
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
        default=False,
        required=False,
        help="Postprocess the output.",
    )
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
        "--dist-threshold",
        type=float,
        required=False,
        help="Distance threshold for generating coordinates.",
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
    # Postprocessing arguments
    parser.add_argument(
        "--postprocess-engine",
        type=str,
        required=False,
        choices=["xtb", "orca"],
        help="QM engine to use for postprocessing.",
    )
    # xTB specific arguments
    parser.add_argument(
        "--xtb-path",
        type=str,
        required=False,
        help="Path to the xTB binary.",
    )
    # ORCA specific arguments
    parser.add_argument(
        "--orca-path",
        type=str,
        required=False,
        help="Path to the ORCA binary.",
    )
    args = parser.parse_args(argv)
    args_dict = vars(args)

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
    }
    rev_args_dict["refine"] = {
        "engine": args_dict["refine_engine"],
        "max_frag_cycles": args_dict["max_frag_cycles"],
    }
    rev_args_dict["generate"] = {
        "min_num_atoms": args_dict["min_num_atoms"],
        "max_num_atoms": args_dict["max_num_atoms"],
        "init_coord_scaling": args_dict["init_coord_scaling"],
        "increase_scaling_factor": args_dict["increase_scaling_factor"],
        "dist_threshold": args_dict["dist_threshold"],
        "element_composition": args_dict["element_composition"],
        "forbidden_elements": args_dict["forbidden_elements"],
    }
    # XTB specific arguments
    rev_args_dict["xtb"] = {"xtb_path": args_dict["xtb_path"]}
    # ORCA specific arguments
    rev_args_dict["orca"] = {"orca_path": args_dict["orca_path"]}
    rev_args_dict["postprocess"] = {"engine": args_dict["postprocess_engine"]}

    return rev_args_dict
