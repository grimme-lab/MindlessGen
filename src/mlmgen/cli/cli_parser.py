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
        "--verbosity",
        type=int,
        required=False,
        choices=[0, 1, 2],
        help="Verbosity level (0, 1, or 2).",
    )
    parser.add_argument(
        "-e",
        "--engine",
        type=str,
        required=False,
        choices=["xtb", "orca"],
        help="QM engine to use.",
    )
    parser.add_argument(
        "-mc",
        "--max-cycles",
        type=int,
        required=False,
        help="Maximum number of optimization cycles.",
    )
    # XTB specific arguments
    # TODO: Add XTB specific arguments
    # ORCA specific arguments
    # TODO: Add ORCA specific arguments
    args = parser.parse_args(argv)
    args_dict = vars(args)

    # General arguments
    rev_args_dict: dict[str, dict] = {}
    rev_args_dict["general"] = {
        "config": args_dict["config"],
        "verbosity": args_dict["verbosity"],
        "engine": args_dict["engine"],
        "max_cycles": args_dict["max_cycles"],
    }
    # XTB specific arguments
    rev_args_dict["xtb"] = {}
    # ORCA specific arguments
    rev_args_dict["orca"] = {}

    return rev_args_dict
