"""
Squarer
=======

Dummy command line tool to square a number.
"""

from .__version__ import __version__
from .cli import console_entry_point
from .generator import generator as square

__all__ = ["__version__", "console_entry_point", "square"]
