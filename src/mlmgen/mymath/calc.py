"""
Mathematical functions.
"""

from __future__ import annotations


def square_a_number(a: float | int) -> float | int:
    if not isinstance(a, (float, int)):
        raise TypeError("Float or int expected.")

    return a * a
