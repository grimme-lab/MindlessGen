from __future__ import annotations

import argparse
from collections.abc import Sequence
from ..__version__ import __version__


def cli_parser(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    """
    # get command line argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=argparse.FileType("r"), help="Input file.", required=False
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args(argv)

    return args
