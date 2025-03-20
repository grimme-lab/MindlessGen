"""
Postprocess the generated molecules.
"""

from threading import Event
from .molecule import Molecule
from ..qm import QMMethod
from ..prog import PostProcessConfig, ResourceMonitor


def postprocess_mol(
    mol: Molecule,
    engine: QMMethod,
    config: PostProcessConfig,
    resources_local: ResourceMonitor,
    stop_event: Event,
    verbosity: int = 1,
) -> Molecule | None:
    """
    Postprocess the generated molecule.

    Arguments:
    mol (Molecule): Molecule to postprocess

    Returns:
    Molecule: Postprocessed molecule
    """

    if config.debug:
        verbosity = 3
    if verbosity > 2:
        print("Postprocessing molecule...")
    if config.optimize:
        try:
            with resources_local.occupy_cores(config.ncores):
                if stop_event.is_set():
                    return None
                postprocmol = engine.optimize(
                    mol,
                    max_cycles=config.opt_cycles,
                    ncores=config.ncores,
                    verbosity=verbosity,
                )
        except RuntimeError as e:
            raise RuntimeError("Optimization in postprocessing failed.") from e
    else:
        try:
            with resources_local.occupy_cores(config.ncores):
                if stop_event.is_set():
                    return None
                engine.singlepoint(mol, config.ncores, verbosity=verbosity)
            postprocmol = mol
        except RuntimeError as e:
            raise RuntimeError(
                f"Single point calculation in postprocessing failed with error: {e}"
            ) from e
    return postprocmol
