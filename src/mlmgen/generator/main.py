"""
Mathematical functions.
"""

from __future__ import annotations
from ..molecules import generate_random_molecule
from ..qm import XTB, get_xtb_path
from ..molecules import postprocess


def generator(config: dict) -> int:
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

    if config["general"]["engine"] == "xtb":
        try:
            xtb_path = get_xtb_path(["xtb_dev", "xtb"])
            if not xtb_path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        engine = XTB(xtb_path, config["general"]["verbosity"])
    else:
        raise NotImplementedError("Engine not implemented.")

    print(f"Config: {config}")
    for cycle in range(config["general"]["max_cycles"]):
        print(f"Cycle {cycle + 1}...")
        #   _____                           _
        #  / ____|                         | |
        # | |  __  ___ _ __   ___ _ __ __ _| |_ ___  _ __
        # | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \| '__|
        # | |__| |  __/ | | |  __/ | | (_| | || (_) | |
        #  \_____|\___|_| |_|\___|_|  \__,_|\__\___/|_|

        mol = generate_random_molecule(config["general"]["verbosity"])

        try:
            #    ____        _   _           _
            #   / __ \      | | (_)         (_)
            #  | |  | |_ __ | |_ _ _ __ ___  _ _______
            #  | |  | | '_ \| __| | '_ ` _ \| |_  / _ \
            #  | |__| | |_) | |_| | | | | | | |/ /  __/
            #   \____/| .__/ \__|_|_| |_| |_|_/___\___|
            #         | |
            #         |_|
            optimized_molecule = postprocess(
                mol=mol, engine=engine, verbosity=config["general"]["verbosity"]
            )
            print("Postprocessing successful. Optimized molecule:")
            print(optimized_molecule)
            optimized_molecule.write_xyz_to_file("optimized_molecule.xyz")
            return 0
        except RuntimeError as e:
            print(f"Postprocessing failed for cycle {cycle + 1}.\n")
            if config["general"]["verbosity"] > 1:
                print(e)
            continue
    raise RuntimeError("Postprocessing failed for all cycles.")
