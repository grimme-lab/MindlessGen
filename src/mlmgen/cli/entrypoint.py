"""
Entrypoint for command line interface.
"""

from __future__ import annotations

from collections.abc import Sequence

from ..generator import generator
from .cli_parser import cli_parser as cl


def console_entry_point(argv: Sequence[str] | None = None) -> int:
    """
    Entrypoint for command line interface.
    """
    args = cl(argv)
    # convert args to dictionary
    kwargs = vars(args)
    raise SystemExit(generator(kwargs))
