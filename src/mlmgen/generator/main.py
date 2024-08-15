"""
Mathematical functions.
"""

from __future__ import annotations

import multiprocessing as mp
import warnings

from ..molecules import generate_random_molecule, Molecule
from ..qm import XTB, get_xtb_path, QMMethod
from ..molecules import postprocess
from ..prog import ConfigManager


def generator(config: ConfigManager) -> Molecule | int | None:
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

    if config.general.print_config:
        print(config)
        return 0

    if config.general.engine == "xtb":
        try:
            xtb_path = get_xtb_path(["xtb_dev", "xtb"])
            if not xtb_path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        engine = XTB(xtb_path, config.general.verbosity)
    else:
        raise NotImplementedError("Engine not implemented.")

    if config.general.verbosity > 0:
        print(config)

    manager = mp.Manager()
    stop_event = manager.Event()
    cycles = range(config.general.max_cycles)

    # lower number of the available cores and the configured parallelism
    num_cores = min(mp.cpu_count(), config.general.parallel)
    if config.general.parallel > mp.cpu_count():
        warnings.warn(
            f"Number of cores requested ({config.general.parallel}) is greater "
            + f"than the number of available cores ({mp.cpu_count()})."
            + f"Using {num_cores} cores instead."
        )

    if config.general.verbosity > 0:
        print(f"Running with {num_cores} cores.")
    backup_verbosity: int | None = None
    if num_cores > 1 and config.general.verbosity > 0:
        backup_verbosity = config.general.verbosity  # Save verbosity level for later
        config.general.verbosity = 0  # Disable verbosity if parallel

    print("Cycle... ", end="", flush=True)
    with mp.Pool(processes=num_cores) as pool:
        results = pool.starmap(
            single_molecule_generator,
            [(config, engine, cycle, stop_event) for cycle in cycles],
        )

    # Restore verbosity level if it was changed
    if backup_verbosity is not None:
        config.general.verbosity = backup_verbosity

    # Filter out None values and return the first successful molecule
    optimized_molecule: Molecule | None = None
    for i, result in enumerate(results):
        if result is not None:
            cycles_needed = i + 1
            optimized_molecule = result

    if optimized_molecule is None:
        raise RuntimeError("Postprocessing failed for all cycles.")

    if config.general.verbosity > 0:
        print(f"\nOptimized mindless molecule found in {cycles_needed} cycles.")
        print(optimized_molecule)
    return optimized_molecule


def single_molecule_generator(
    config: ConfigManager, engine: QMMethod, cycle: int, stop_event
):
    """
    Generate a single molecule.
    """
    if stop_event.is_set():
        return None  # Exit early if a molecule has already been found

    if config.general.verbosity == 0:
        # print the cycle in one line, not starting a new line
        print("✔", end="", flush=True)
    else:
        print(f"Cycle {cycle + 1}:")
    #   _____                           _
    #  / ____|                         | |
    # | |  __  ___ _ __   ___ _ __ __ _| |_ ___  _ __
    # | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \| '__|
    # | |__| |  __/ | | |  __/ | | (_| | || (_) | |
    #  \_____|\___|_| |_|\___|_|  \__,_|\__\___/|_|

    mol = generate_random_molecule(config.general)

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
            mol=mol,
            engine=engine,
            min_nat=config.general.min_num_atoms,
            max_nat=config.general.max_num_atoms,
            verbosity=config.general.verbosity,
        )
        if not stop_event.is_set():
            stop_event.set()  # Signal other processes to stop
            if config.general.verbosity > 1:
                print("Postprocessing successful. Optimized molecule:")
                print(optimized_molecule)
            return optimized_molecule
    except RuntimeError as e:
        if config.general.verbosity > 0:
            print(f"Postprocessing failed for cycle {cycle + 1}.")
            if config.general.verbosity > 1:
                print(e)
        return None
