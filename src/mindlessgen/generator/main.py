"""
Main driver of MindlessGen.
"""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import Future, as_completed, wait
from pathlib import Path
import multiprocessing as mp
from threading import Event
import warnings
from dataclasses import dataclass

from ..molecules import generate_random_molecule, Molecule
from ..qm import XTB, get_xtb_path, QMMethod, ORCA, get_orca_path, GXTB, get_gxtb_path
from ..molecules import iterative_optimization, postprocess_mol
from ..prog import ConfigManager, ParallelManager
from ..prog.config import MINCORES_PLACEHOLDER

from ..__version__ import __version__

MINDLESS_MOLECULES_FILE = "mindless.molecules"


@dataclass
class Block:
    num_molecules: int
    ncores: int


def generator(config: ConfigManager) -> tuple[list[Molecule], int]:
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

    if config.general.verbosity > 0:
        print(header(str(__version__)))

    if config.general.print_config:
        print(config)
        return [], 0

    # Import and set up required engines
    refine_engine: QMMethod = setup_engines(
        config.refine.engine,
        config,
        get_xtb_path,
        get_orca_path,  # g-xTB cannot be used anyway
    )

    if config.general.postprocess:
        postprocess_engine: QMMethod | None = setup_engines(
            config.postprocess.engine,
            config,
            get_xtb_path,
            get_orca_path,
            get_gxtb_path,
        )
    else:
        postprocess_engine = None

    if config.general.verbosity > 0:
        print(config)

    config.check_config(verbosity=config.general.verbosity)

    num_cores = min(mp.cpu_count(), config.general.parallel)
    if config.general.verbosity > 0:
        print(f"Running with {num_cores} cores.")

    # Check if the file "mindless.molecules" exists. If yes, append to it.
    if Path(MINDLESS_MOLECULES_FILE).is_file():
        if config.general.verbosity > 0:
            print(f"\n--- Appending to existing file '{MINDLESS_MOLECULES_FILE}'. ---")

    backup_verbosity: int | None = None
    if num_cores > 1 and config.general.verbosity > 0:
        backup_verbosity = config.general.verbosity  # Save verbosity level for later
        config.general.verbosity = 0  # Disable verbosity if parallel

    exitcode = 0
    optimized_molecules: list[Molecule] = []

    # Initialize parallel blocks here
    blocks = setup_blocks(num_cores, config.general.num_molecules)
    blocks.sort(key=lambda x: x.ncores)

    # Set up parallel blocks environment
    with ParallelManager(num_cores // MINCORES_PLACEHOLDER, num_cores) as parallel:
        # The following creates a queue of futures which occupy a certain number of cores each
        # as defined by each block
        # Each future represents the generation of one molecule
        # NOTE: proceeding this way assures that each molecule gets a static number of cores
        # a dynamic setting would also be thinkable and straightforward to implement
        tasks: list[Future[Molecule | None]] = []
        for block in blocks:
            for _ in range(block.num_molecules):
                # TODO: remove objects from future call that cannot be pickled (basically everything that doesn't consist only of primitives)
                tasks.append(
                    parallel.executor.submit(
                        single_molecule_generator,
                        len(tasks) + 1,
                        config,
                        parallel,
                        refine_engine,
                        postprocess_engine,
                        block.ncores,
                    )
                )

        # Collect results of all tries to create a molecule
        results: list[Molecule | None] = [task.result() for task in as_completed(tasks)]

        # Restore verbosity level if it was changed
        if backup_verbosity is not None:
            config.general.verbosity = backup_verbosity

        for molcount, optimized_molecule in enumerate(results):
            if optimized_molecule is None:
                # TODO: molcount might not align with the number of the molecule that actually failed, look into this
                warnings.warn(
                    "Molecule generation including optimization (and postprocessing) "
                    + f"failed for all cycles for molecule {molcount + 1}."
                )
                exitcode = 1
                continue

            # if config.general.verbosity > 0:
            #     print(f"Optimized mindless molecule found in {cycles_needed} cycles.")
            #     print(optimized_molecule)

            if config.general.write_xyz:
                optimized_molecule.write_xyz_to_file()
                if config.general.verbosity > 0:
                    print(
                        f"Written molecule file 'mlm_{optimized_molecule.name}.xyz'.\n"
                    )
                with open("mindless.molecules", "a", encoding="utf8") as f:
                    f.write(f"mlm_{optimized_molecule.name}\n")

            optimized_molecules.append(optimized_molecule)

    return optimized_molecules, exitcode


def single_molecule_generator(
    molcount: int,
    config: ConfigManager,
    parallel: ParallelManager,
    refine_engine: QMMethod,
    postprocess_engine: QMMethod | None,
    ncores: int,
) -> Molecule | None:
    """
    Generate a single molecule (from start to finish).
    """

    # Wait for enough cores
    parallel.occupy_cores(ncores)

    # print a decent header for each molecule iteration
    # if config.general.verbosity > 0:
    #     print(f"\n{'='*80}")
    #     print(
    #         f"{'='*22} Generating molecule {molcount + 1:<4} of "
    #         + f"{config.general.num_molecules:<4} {'='*24}"
    #     )
    #     print(f"{'='*80}")

    with ParallelManager(ncores, ncores) as parallel_local:
        stop_event = parallel_local.manager.Event()

        # Launch worker processes to find molecule
        # if config.general.verbosity == 0:
        #     print("Cycle... ", end="", flush=True)
        cycles = range(config.general.max_cycles)
        tasks: list[Future[Molecule | None]] = []
        for cycle in cycles:
            # TODO: remove objects from future call that cannot be pickled (basically everything that doesn't consist only of primitives)
            tasks.append(
                parallel_local.executor.submit(
                    single_molecule_step,
                    config,
                    parallel_local,
                    refine_engine,
                    postprocess_engine,
                    cycle,
                    stop_event,
                )
            )

        # Finally, add a future to set the stop_event if all jobs are completed
        # parallel_local.executor.submit(
        #     lambda: stop_event.set() if wait(tasks) else None
        # )
        #
        # stop_event.wait()

        results = [task.result() for task in as_completed(tasks)]

    # if config.general.verbosity == 0:
    #     print("")

    optimized_molecule: Molecule | None = None
    for i, result in enumerate(results):
        if result is not None:
            cycles_needed = i + 1
            optimized_molecule = result
            break

    # if config.general.verbosity > 0:
    #     print(f"Optimized mindless molecule found in {cycles_needed} cycles.")
    #     print(optimized_molecule)

    # Free up the cores
    parallel.free_cores(ncores)

    return optimized_molecule


def single_molecule_step(
    config: ConfigManager,
    parallel: ParallelManager,
    refine_engine: QMMethod,
    postprocess_engine: QMMethod | None,
    cycle: int,
    stop_event: Event,
) -> Molecule | None:
    """Execute one step in a single molecule generation"""

    if stop_event.is_set():
        return None  # Exit early if a molecule has already been found

    # if config.general.verbosity == 0:
    #     # print the cycle in one line, not starting a new line
    #     print("✔", end="", flush=True)
    # elif config.general.verbosity > 0:
    #     print(f"Cycle {cycle + 1}:")

    #   _____                           _
    #  / ____|                         | |
    # | |  __  ___ _ __   ___ _ __ __ _| |_ ___  _ __
    # | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \| '__|
    # | |__| |  __/ | | |  __/ | | (_| | || (_) | |
    #  \_____|\___|_| |_|\___|_|  \__,_|\__\___/|_|

    try:
        mol = generate_random_molecule(config.generate, config.general.verbosity)
    # RuntimeError is not caught here, as in this part, runtime errors are not expected to occur
    # and shall therefore be raised to the main function
    except (
        SystemExit
    ) as e:  # debug functionality: raise SystemExit to stop the whole execution
        # if config.general.verbosity > 0:
        #     print(f"Generation aborted for cycle {cycle + 1}.")
        #     if config.general.verbosity > 1:
        #         print(e)
        stop_event.set()
        return None
    except RuntimeError as e:
        # if config.general.verbosity > 0:
        #     print(f"Generation failed for cycle {cycle + 1}.")
        #     if config.general.verbosity > 1:
        #         print(e)
        return None

    try:
        #    ____        _   _           _
        #   / __ \      | | (_)         (_)
        #  | |  | |_ __ | |_ _ _ __ ___  _ _______
        #  | |  | | '_ \| __| | '_ ` _ \| |_  / _ \
        #  | |__| | |_) | |_| | | | | | | |/ /  __/
        #   \____/| .__/ \__|_|_| |_| |_|_/___\___|
        #         | |
        #         |_|
        optimized_molecule = iterative_optimization(
            mol,
            refine_engine,
            config.generate,
            config.refine,
            parallel,
            verbosity=config.general.verbosity,
        )
    except RuntimeError as e:
        # if config.general.verbosity > 0:
        #     print(f"Refinement failed for cycle {cycle + 1}.")
        #     if config.general.verbosity > 1 or config.refine.debug:
        #         print(e)
        return None
    finally:
        if config.refine.debug:
            stop_event.set()

    if config.general.postprocess:
        try:
            optimized_molecule = postprocess_mol(
                optimized_molecule,
                postprocess_engine,  # type: ignore
                config.postprocess,
                parallel,
                verbosity=config.general.verbosity,
            )
        except RuntimeError as e:
            # if config.general.verbosity > 0:
            #     print(f"Postprocessing failed for cycle {cycle + 1}.")
            #     if config.general.verbosity > 1 or config.postprocess.debug:
            #         print(e)
            return None
        finally:
            if config.postprocess.debug:
                stop_event.set()  # Stop further runs if debugging of this step is enabled
        # if config.general.verbosity > 1:
        #     print("Postprocessing successful.")

    if not stop_event.is_set():
        stop_event.set()  # Signal other processes to stop
        return optimized_molecule
    elif config.refine.debug or config.postprocess.debug:
        return optimized_molecule
    else:
        return None


def header(version: str) -> str:
    """
    This function prints the header of the program.
    """
    headerstr = (
        # pylint: disable=C0301
        "╔══════════════════════════════════════════════════════════════════════════════════════════════════╗\n"
        "║                                                                                                  ║\n"
        "║   ███╗   ███╗██╗███╗   ██╗██████╗ ██╗     ███████╗███████╗███████╗ ██████╗ ███████╗███╗   ██╗    ║\n"
        "║   ████╗ ████║██║████╗  ██║██╔══██╗██║     ██╔════╝██╔════╝██╔════╝██╔════╝ ██╔════╝████╗  ██║    ║\n"
        "║   ██╔████╔██║██║██╔██╗ ██║██║  ██║██║     █████╗  ███████╗███████╗██║  ███╗█████╗  ██╔██╗ ██║    ║\n"
        "║   ██║╚██╔╝██║██║██║╚██╗██║██║  ██║██║     ██╔══╝  ╚════██║╚════██║██║   ██║██╔══╝  ██║╚██╗██║    ║\n"
        "║   ██║ ╚═╝ ██║██║██║ ╚████║██████╔╝███████╗███████╗███████║███████║╚██████╔╝███████╗██║ ╚████║    ║\n"
        "║   ╚═╝     ╚═╝╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚══════╝╚══════╝╚══════╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝    ║\n"
        "║                                                                                                  ║\n"
        f"║                                       MindlessGen v{version[:5]}                                         ║\n"
        "║                                 Semi-Automated Molecule Generator                                ║\n"
        "║                                                                                                  ║\n"
        "║                          Licensed under the Apache License, Version 2.0                          ║\n"
        "║                           (http://www.apache.org/licenses/LICENSE-2.0)                           ║\n"
        "╚══════════════════════════════════════════════════════════════════════════════════════════════════╝"
    )
    return headerstr


# Define a utility function to set up the required engine
def setup_engines(
    engine_type: str,
    cfg: ConfigManager,
    xtb_path_func: Callable,
    orca_path_func: Callable,
    gxtb_path_func: Callable | None = None,
):
    """
    Set up the required engine.
    """
    if engine_type == "xtb":
        try:
            path = xtb_path_func(cfg.xtb.xtb_path)
            if not path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        return XTB(path, cfg.xtb)
    elif engine_type == "orca":
        try:
            path = orca_path_func(cfg.orca.orca_path)
            if not path:
                raise ImportError("orca not found.")
        except ImportError as e:
            raise ImportError("orca not found.") from e
        return ORCA(path, cfg.orca)
    elif engine_type == "gxtb":
        if gxtb_path_func is None:
            raise ImportError("No callable function for determining the g-xTB path.")
        path = gxtb_path_func(cfg.gxtb.gxtb_path)
        if not path:
            raise ImportError("'gxtb' binary could not be found.")
        return GXTB(path, cfg.gxtb)
    else:
        raise NotImplementedError("Engine not implemented.")


def setup_blocks(ncores: int, num_molecules: int) -> list[Block]:
    blocks: list[Block] = []

    # Maximum and minimum number of parallel processes possible
    maxcores = ncores
    mincores = MINCORES_PLACEHOLDER
    maxprocs = max(1, ncores // mincores)
    minprocs = max(1, ncores // maxcores)

    # Distribute number of molecules among blocks
    # First (if possible) create the maximum number of parallel blocks (maxprocs) and distribute as many molecules as possible
    molecules_left = num_molecules
    if molecules_left >= maxprocs:
        p = maxprocs
        molecules_per_block = molecules_left // p
        for _ in range(p):
            blocks.append(Block(molecules_per_block, ncores // p))
        molecules_left -= molecules_per_block * p

    # While there are more than minprocs (1) molecules left find the optimal number of parallel blocks
    # Again distribute as many molecules per block as possible
    while molecules_left >= minprocs:
        p = max(
            [
                j
                for j in range(minprocs, maxprocs)
                if ncores % j == 0 and j <= molecules_left
            ]
        )
        molecules_per_block = molecules_left // p
        for _ in range(p):
            blocks.append(Block(molecules_per_block, ncores // p))
        molecules_left -= molecules_per_block * p

    # NOTE: using minprocs = 1 this is probably never true
    if molecules_left > 0:
        blocks.append(Block(molecules_left, maxcores))
        molecules_left -= molecules_left

    return blocks
