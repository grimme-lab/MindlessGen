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
    bool: Success status
    """

    # Run a singlepoint calculation with the postprocess engine to see if it is stable
    try:
        engine.singlepoint(mol, verbosity=verbosity)
    except RuntimeError as e:
        raise RuntimeError("Single point calculation in postprocessing failed.") from e
    return mol
