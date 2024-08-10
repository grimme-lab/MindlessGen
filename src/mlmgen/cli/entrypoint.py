"""
Entrypoint for command line interface.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from ..mymath import square_a_number as square


def console_entry_point(argv: Sequence[str] | None = None) -> int:
    # get command line argument
    parser = argparse.ArgumentParser()
    parser.add_argument("number", type=float, help="Number to square.")
    args = parser.parse_args(argv)

    # print result
    print(square(args.number))

    return 0
