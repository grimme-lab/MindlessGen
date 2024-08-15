"""
Command line interface
======================

This module contains functionality for the CLI.
"""

from .entrypoint import console_entry_point
from .cli_parser import cli_parser

__all__ = ["console_entry_point", "cli_parser"]
