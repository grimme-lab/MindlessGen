"""
Test the xtb interface.
"""

from __future__ import annotations
from pathlib import Path
import pytest
import numpy as np
from mindlessgen.qm import XTB, get_xtb_path  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore
from mindlessgen.prog import XTBConfig  # type: ignore


# mark all tests as optional as they depend on the availability of xtb
# test the XTB optimizer using the ethanol molecule
@pytest.mark.optional
def test_xtb_optimize_xtb(coordinates_ethanol: np.ndarray) -> None:
    """
    Test the optimization of ethanol with xtb.
    """
    cfg = XTBConfig()
    xtb_path = get_xtb_path()
    if xtb_path:
        xtb = XTB(xtb_path, cfg)
    else:
        raise RuntimeError("xtb not found.")
    mol = Molecule()
    mol.xyz = coordinates_ethanol
    mol.ati = np.array([7, 5, 5, 0, 0, 0, 0, 0, 0])
    mol.charge = 0
    mol.uhf = 0
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
        atol=1e-4,
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


# load the molecule: C2H4N1O1Au1
@pytest.fixture
def mol_C2H4N1O1Au1() -> Molecule:
    """
    Load the molecule C2H4N1O1Au1 from fixtues/C2H4N1O1Au1.xyz.
    """
    mol = Molecule("C2H4N1O1Au1")
    # get the Path of this file
    testsdir = Path(__file__).resolve().parents[1]
    xyz_file = testsdir / "fixtures/C2H4N1O1Au1.xyz"
    mol.read_xyz_from_file(xyz_file)
    with open(testsdir / "fixtures/C2H4N1O1Au1.CHRG", encoding="utf8") as f:
        mol.charge = int(f.read())
    mol.uhf = 0
    return mol


# load the molecule: H3B4Pd1Rn1
@pytest.fixture
def mol_H3B4Pd1Rn1() -> Molecule:
    """
    Load the molecule H3B4Pd1Rn1 from fixtues/H3B4Pd1Rn1.xyz.
    """
    mol = Molecule("H3B4Pd1Rn1")
    # get the Path of this file
    testsdir = Path(__file__).resolve().parents[1]
    xyz_file = testsdir / "fixtures/H3B4Pd1Rn1.xyz"
    mol.read_xyz_from_file(xyz_file)
    with open(testsdir / "fixtures/H3B4Pd1Rn1.CHRG", encoding="utf8") as f:
        mol.charge = int(f.read())
    mol.uhf = 0
    return mol


@pytest.mark.optional
def test_check_gap_low_gap(mol_C2H4N1O1Au1: Molecule, mol_H3B4Pd1Rn1: Molecule):
    """
    Test check_gap with a molecule that has a low HOMO-LUMO gap and one with a sufficient gap.
    """

    cfg = XTBConfig()
    try:
        xtb_path = get_xtb_path()
        if not xtb_path:
            raise ImportError("xtb not found.")
    except ImportError as e:
        raise ImportError("xtb not found.") from e
    engine = XTB(xtb_path, cfg)
    # Test for molecule with low gap
    result_low_gap = engine.check_gap(mol_H3B4Pd1Rn1, threshold=0.5)
    assert result_low_gap is False

    # Test for molecule with high gap
    result_high_gap = engine.check_gap(mol_C2H4N1O1Au1, threshold=0.5)
    assert result_high_gap is True
