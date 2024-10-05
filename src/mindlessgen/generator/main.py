"""
Main driver of MindlessGen.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
import multiprocessing as mp
import warnings

from ..molecules import generate_random_molecule, Molecule
from ..qm import XTB, get_xtb_path, QMMethod, ORCA, get_orca_path, GP3, get_gp3_path
from ..molecules import iterative_optimization, postprocess_mol
from ..prog import ConfigManager

from .. import __version__

MINDLESS_MOLECULES_FILE = "mindless.molecules"


def generator(config: ConfigManager) -> tuple[list[Molecule] | None, int]:
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
        return None, 0

    # Import and set up required engines
    refine_engine: QMMethod = setup_engines(
        config.refine.engine,
        config,
        get_xtb_path,
        get_orca_path,  # GP3 cannot be used anyway
    )

    if config.general.postprocess:
        postprocess_engine: QMMethod | None = setup_engines(
            config.postprocess.engine, config, get_xtb_path, get_orca_path, get_gp3_path
        )
    else:
        postprocess_engine = None

    if config.general.verbosity > 0:
        print(config)

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

    if num_cores > 1 and config.general.verbosity > 0:
        # raise warning that parallelization will disable verbosity
        warnings.warn(
            "Parallelization will disable verbosity during iterative search. "
            + "Set '--verbosity 0' or '-P 1' to avoid this warning, or simply ignore it."
        )
    if num_cores > 1 and config.postprocess.debug:
        # raise warning that debugging of postprocessing will disable parallelization
        warnings.warn(
            "Debug output might seem to be redundant due to the parallel processes "
            + "with possibly similar errors in parallel mode. "
            + "Don't be confused!"
        )

    # Check if the file "mindless.molecules" exists. If yes, append to it.
    if Path(MINDLESS_MOLECULES_FILE).is_file():
        if config.general.verbosity > 0:
            print(f"\n--- Appending to existing file '{MINDLESS_MOLECULES_FILE}'. ---")
    exitcode = 0
    optimized_molecules: list[Molecule] = []
    for molcount in range(config.general.num_molecules):
        # print a decent header for each molecule iteration
        if config.general.verbosity > 0:
            print(f"\n{'='*80}")
            print(
                f"{'='*22} Generating molecule {molcount + 1:<4} of "
                + f"{config.general.num_molecules:<4} {'='*24}"
            )
            print(f"{'='*80}")
        manager = mp.Manager()
        stop_event = manager.Event()
        cycles = range(config.general.max_cycles)
        backup_verbosity: int | None = None
        if num_cores > 1 and config.general.verbosity > 0:
            backup_verbosity = (
                config.general.verbosity
            )  # Save verbosity level for later
            config.general.verbosity = 0  # Disable verbosity if parallel

        if config.general.verbosity == 0:
            print("Cycle... ", end="", flush=True)
        with mp.Pool(processes=num_cores) as pool:
            results = pool.starmap(
                single_molecule_generator,
                [
                    (config, refine_engine, postprocess_engine, cycle, stop_event)
                    for cycle in cycles
                ],
            )
        if config.general.verbosity == 0:
            print("")

        # Restore verbosity level if it was changed
        if backup_verbosity is not None:
            config.general.verbosity = backup_verbosity

        # Filter out None values and return the first successful molecule
        optimized_molecule: Molecule | None = None
        for i, result in enumerate(results):
            if result is not None:
                cycles_needed = i + 1
                optimized_molecule = result
                break

        if optimized_molecule is None:
            warnings.warn(
                "Molecule generation including optimization (and postprocessing) "
                + f"failed for all cycles for molecule {molcount + 1}."
            )
            exitcode = 1
            continue
        if config.general.verbosity > 0:
            print(f"Optimized mindless molecule found in {cycles_needed} cycles.")
            print(optimized_molecule)
        if config.general.write_xyz:
            optimized_molecule.write_xyz_to_file()
            if config.general.verbosity > 0:
                print(f"Written molecule file 'mlm_{optimized_molecule.name}.xyz'.\n")
            with open("mindless.molecules", "a", encoding="utf8") as f:
                f.write(f"mlm_{optimized_molecule.name}\n")
        optimized_molecules.append(optimized_molecule)

    return optimized_molecules, exitcode


def single_molecule_generator(
    config: ConfigManager,
    refine_engine: QMMethod,
    postprocess_engine: QMMethod | None,
    cycle: int,
    stop_event,
) -> Molecule | None:
    """
    Generate a single molecule.
    """
    if stop_event.is_set():
        return None  # Exit early if a molecule has already been found

    if config.general.verbosity == 0:
        # print the cycle in one line, not starting a new line
        print("✔", end="", flush=True)
    elif config.general.verbosity > 0:
        print(f"Cycle {cycle + 1}:")
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
        if config.general.verbosity > 0:
            print(f"Generation aborted for cycle {cycle + 1}.")
            if config.general.verbosity > 1:
                print(e)
        stop_event.set()
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
            mol=mol,
            engine=refine_engine,
            config_generate=config.generate,
            config_refine=config.refine,
            verbosity=config.general.verbosity,
        )
    except RuntimeError as e:
        if config.general.verbosity > 0:
            print(f"Refinement failed for cycle {cycle + 1}.")
            if config.general.verbosity > 1 or config.refine.debug:
                print(e)
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
                config.general.verbosity,
            )
        except RuntimeError as e:
            if config.general.verbosity > 0:
                print(f"Postprocessing failed for cycle {cycle + 1}.")
                if config.general.verbosity > 1 or config.postprocess.debug:
                    print(e)
            return None
        finally:
            if config.postprocess.debug:
                stop_event.set()  # Stop further runs if debugging of this step is enabled
        if config.general.verbosity > 1:
            print("Postprocessing successful.")

    if config.general.gxtb_development:
        # additional g-xTB engine after refinement and postprocessing for development purposes
        gp3 = GP3(get_gp3_path())
        try:
            _gxtb_dev_check(
                optimized_molecule,
                gp3,
                config.general.gxtb_scf_cycles,
                config.general.verbosity,
            )
        except (RuntimeError, ValueError) as e:
            if config.general.verbosity > 0:
                print(f"g-xTB postprocessing failed for cycle {cycle + 1}.")
                if config.general.verbosity > 1:
                    print(e)
            return None
        if config.general.verbosity > 1:
            print("g-xTB postprocessing successful.")

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
    gp3_path_func: Callable | None = None,
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
    elif engine_type == "gp3":
        if gp3_path_func is None:
            raise ImportError("No callable function for determining the gp3 path.")
        path = gp3_path_func()
        if not path:
            raise ImportError("'gp3' binary could not be found.")
        return GP3(path)
    else:
        raise NotImplementedError("Engine not implemented.")


def _gxtb_dev_check(
    mol: Molecule, gp3: GP3, scf_iter_limit: int, verbosity: int = 0
) -> None:
    """
    ONLY FOR IN-HOUSE g-xTB DEVELOPMENT PURPOSES: Check the SCF iterations of the cation and anion.
    """
    # 1) Single point calculation with g-xTB for the cation
    tmp_mol = mol.copy()
    tmp_mol.charge += 1
    tmp_mol.uhf += 1
    if mol.uhf != 0 and verbosity > 0:
        print(
            f"WARNING: UHF value was set to {mol.uhf} and not 0. "
            + "For the g-xTB cationic calculation, we are increasing it by 1. "
            + "(Could be ill-defined.)"
        )
    gxtb_output = gp3.singlepoint(tmp_mol, verbosity=verbosity)
    # gp3_output looks like this:
    # [...]
    #   13     -155.03101038        0.00000000        0.00000001       16.45392733   8    F
    #           13  scf iterations
    #           eigenvalues
    # [...]
    # Check for the number of scf iterations
    scf_iterations = 0
    for line in gxtb_output.split("\n"):
        if "scf iterations" in line:
            scf_iterations = int(line.split()[0])
            break
    if scf_iterations == 0:
        raise ValueError("SCF iterations not found in GP3 output.")
    if scf_iterations > scf_iter_limit:
        raise ValueError(f"SCF iterations exceeded limit of {scf_iter_limit}.")

    # 2) Single point calculation with g-xTB for the anion
    tmp_mol = mol.copy()
    tmp_mol.charge -= 1
    tmp_mol.uhf += 1
    if mol.uhf != 0 and verbosity > 0:
        print(
            f"WARNING: UHF value was set to {mol.uhf} and not 0. "
            + "For the g-xTB anionic calculation, we are increasing it by 1. "
            + "(Could be ill-defined.)"
        )
    gxtb_output = gp3.singlepoint(tmp_mol, verbosity=verbosity)
    # Check for the number of scf iterations
    scf_iterations = 0
    for line in gxtb_output.split("\n"):
        if "scf iterations" in line:
            scf_iterations = int(line.split()[0])
            break
    if scf_iterations == 0:
        raise ValueError("SCF iterations not found in GP3 output.")
    if scf_iterations > scf_iter_limit:
        raise ValueError(f"SCF iterations exceeded limit of {scf_iter_limit}.")
