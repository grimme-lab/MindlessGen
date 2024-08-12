"""
Mathematical functions.
"""

from __future__ import annotations


from ..molecules import generate_molecule


def generator(input: str | None = None, verbosity: int = 1) -> int:
    """
    Generate a molecule.
    """
    if input:
        print(f"Input file: {input}")
    else:
        mol = generate_molecule(verbosity)
        mol.write_xyz_to_file("molecule.xyz")

    return 0
