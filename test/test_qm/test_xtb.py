"""
Test the xtb interface.
"""

from __future__ import annotations
import pytest
import numpy as np
from mlmgen.qm import XTB, get_xtb_path  # type: ignore
from mlmgen.molecules import Molecule  # type: ignore


# mark all tests as optional as they depend on the availability of xtb
# test the XTB optimizer using the ethanol molecule
@pytest.mark.optional
def test_xtb_optimize_xtb(coordinates_ethanol: np.ndarray) -> None:
    """
    Test the optimization of ethanol with xtb.
    """
    xtb_path = get_xtb_path(["xtb_dev", "xtb"])
    if xtb_path:
        xtb = XTB(xtb_path)
    else:
        raise RuntimeError("xtb not found.")
    mol = Molecule()
    mol.xyz = coordinates_ethanol
    mol.ati = np.array([7, 5, 5, 0, 0, 0, 0, 0, 0])
    mol.charge = 0
    mol.num_atoms = 9

    optimized_molecule = xtb.optimize(mol)

    print(optimized_molecule)

    # assert that the optimized molecule corresponds approximately to the
    # following coordinates:
    assert np.allclose(
        optimized_molecule.xyz,
        np.array(
            [
                [-1.20859513212472, 0.25336962100001, 0.00620969500040],
                [-0.04956307030717, -0.54667263216476, 0.01643723111179],
                [1.23544832549706, 0.27470686717811, 0.00514725943174],
                [-0.05587883585542, -1.21527863145539, 0.89070772088766],
                [-0.11511912857037, -1.15968668490666, -0.88527523735804],
                [2.09784948781593, -0.38451363822015, -0.04069118001598],
                [1.24929017417545, 0.93458059063697, -0.85902012085995],
                [1.30867409034452, 0.88124403031176, 0.90705730478035],
                [-1.17300591097528, 0.85675047762011, 0.75712732702203],
            ]
        ),
        rtol=1e-3,
    )


@pytest.fixture
def coordinates_ethanol() -> np.ndarray:
    """
    Return the coordinates of ethanol.
    """
    return np.array(
        [
            [-1.1712, 0.2997, 0.0],
            [-0.0463, -0.5665, 0.0],
            [1.2175, 0.2668, 0.0],
            [-0.0958, -1.212, 0.8819],
            [-0.0952, -1.1938, -0.8946],
            [2.105, -0.372, -0.0177],
            [1.2426, 0.9307, -0.8704],
            [1.2616, 0.9052, 0.8886],
            [-1.1291, 0.8364, 0.8099],
        ],
        dtype=float,
    )
