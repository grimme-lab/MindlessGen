"""
Mathematical functions.
"""

from __future__ import annotations


from ..molecules import generate_molecule


def generator(input: str | None = None, verbosity: int | None = 1) -> int:
    """
    Generate a molecule.
    """
    if input is not None:
        print(f"Input file: {input}")
    else:
        natoms = generate_molecule(verbosity)
    print(f"Generated molecule:\n{natoms}")
    return 0
