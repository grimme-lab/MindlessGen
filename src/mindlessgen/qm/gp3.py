"""
This module contains all interactions with the GP3-xTB binary for next-gen tight-binding calculations.
"""

from .base import QMMethod


class GP3(QMMethod):
    """
    This class handles all interaction with the GP3 external dependency.
    """
