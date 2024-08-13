"""
Mathematical functions.
"""

from __future__ import annotations


from ..molecules import generate_molecule
from ..qm import XTB, get_xtb_path


def generator(input: str | None = None, verbosity: int = 1) -> int:
    """
    Generate a molecule.
    """
    xtb_path = get_xtb_path(["xtb_dev", "xtb"])
    if xtb_path:
        xtb = XTB(xtb_path, verbosity)
    else:
        raise RuntimeError("xtb not found.")

    if input:
        print(f"Input file: {input}")
    else:
        mol = generate_molecule(verbosity)

    optimized_molecule = xtb.optimize(mol)
    optimized_molecule.write_xyz_to_file("optimized_molecule.xyz")

    return 0
