"""
Mathematical functions.
"""

from __future__ import annotations
from ..molecules import generate_random_molecule
from ..qm import XTB, get_xtb_path


def generator(inputdict: dict) -> int:
    """
    Generate a molecule.
    """
    # Get desired engine from inputdict
    if inputdict["engine"] == "xtb":
        try:
            xtb_path = get_xtb_path(["xtb_dev", "xtb"])
            if not xtb_path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        xtb = XTB(xtb_path, inputdict["verbosity"])
    else:
        raise NotImplementedError("Engine not implemented.")

    if inputdict["input"]:
        print(f"Input file: {input}")
    else:
        mol = generate_random_molecule(inputdict["verbosity"])

    optimized_molecule = xtb.optimize(mol)
    optimized_molecule.write_xyz_to_file("optimized_molecule.xyz")

    return 0
