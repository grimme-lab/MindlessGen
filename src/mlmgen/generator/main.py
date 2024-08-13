"""
Mathematical functions.
"""

from __future__ import annotations
from ..molecules import generate_random_molecule
from ..qm import XTB, get_xtb_path
from ..molecules import postprocess


def generator(inputdict: dict) -> int:
    """
    Generate a molecule.
    """

    #  ______             _
    # |  ____|           (_)
    # | |__   _ __   __ _ _ _ __   ___
    # |  __| | '_ \ / _` | | '_ \ / _ \
    # | |____| | | | (_| | | | | |  __/)
    # |______|_| |_|\__, |_|_| |_|\___|
    #                __/ |
    #               |___/

    if inputdict["engine"] == "xtb":
        try:
            xtb_path = get_xtb_path(["xtb_dev", "xtb"])
            if not xtb_path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        engine = XTB(xtb_path, inputdict["verbosity"])
    else:
        raise NotImplementedError("Engine not implemented.")

    #   _____                           _
    #  / ____|                         | |
    # | |  __  ___ _ __   ___ _ __ __ _| |_ ___  _ __
    # | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \| '__|
    # | |__| |  __/ | | |  __/ | | (_| | || (_) | |
    #  \_____|\___|_| |_|\___|_|  \__,_|\__\___/|_|

    if inputdict["input"]:
        print(f"Input file: {input}")
    else:
        mol = generate_random_molecule(inputdict["verbosity"])

    #    ____        _   _           _
    #   / __ \      | | (_)         (_)
    #  | |  | |_ __ | |_ _ _ __ ___  _ _______
    #  | |  | | '_ \| __| | '_ ` _ \| |_  / _ \
    #  | |__| | |_) | |_| | | | | | | |/ /  __/
    #   \____/| .__/ \__|_|_| |_| |_|_/___\___|
    #         | |
    #         |_|
    optimized_molecule = postprocess(
        mol=mol, engine=engine, verbosity=inputdict["verbosity"]
    )
    optimized_molecule.write_xyz_to_file("optimized_molecule.xyz")
    return 0
