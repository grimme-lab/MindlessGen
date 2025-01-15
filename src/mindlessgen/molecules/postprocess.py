"""
Postprocess the generated molecules.
"""

from .molecule import Molecule
from ..qm import QMMethod
from ..prog import PostProcessConfig, ParallelManager
from ..prog.config import MINCORES_PLACEHOLDER


def postprocess_mol(
    mol: Molecule,
    engine: QMMethod,
    config: PostProcessConfig,
    parallel: ParallelManager,
    verbosity: int = 1,
) -> Molecule:
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
            parallel.occupy_cores(MINCORES_PLACEHOLDER)
            postprocmol = engine.optimize(
                mol, max_cycles=config.opt_cycles, verbosity=verbosity
            )
        except RuntimeError as e:
            raise RuntimeError("Optimization in postprocessing failed.") from e
        finally:
            parallel.free_cores(MINCORES_PLACEHOLDER)
    else:
        try:
            parallel.occupy_cores(MINCORES_PLACEHOLDER)
            engine.singlepoint(mol, verbosity=verbosity)
            postprocmol = mol
        except RuntimeError as e:
            raise RuntimeError(
                "Single point calculation in postprocessing failed."
            ) from e
        finally:
            parallel.free_cores(MINCORES_PLACEHOLDER)
    return postprocmol
