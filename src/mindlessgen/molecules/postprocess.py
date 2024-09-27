"""
Postprocess the generated molecules.
"""

from .molecule import Molecule
from ..qm import QMMethod
from ..prog import PostProcessConfig


def postprocess_mol(
    mol: Molecule, engine: QMMethod, config: PostProcessConfig, verbosity: int = 1
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
            postprocmol = engine.optimize(
                mol, max_cycles=config.opt_cycles, verbosity=verbosity
            )
        except RuntimeError as e:
            raise RuntimeError("Optimization in postprocessing failed.") from e
    else:
        try:
            engine.singlepoint(mol, verbosity=verbosity)
            postprocmol = mol
        except RuntimeError as e:
            raise RuntimeError(
                "Single point calculation in postprocessing failed."
            ) from e
    return postprocmol
