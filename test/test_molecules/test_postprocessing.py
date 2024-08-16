"""
Test the postprocessing functions in the molecules module.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest
from mindlessgen.prog import ConfigManager  # type: ignore
from mindlessgen.molecules import detect_fragments  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore
from mindlessgen.molecules import iterative_optimization  # type: ignore
from mindlessgen.qm import XTB, get_xtb_path  # type: ignore


@pytest.fixture
def mol_H6O2B2Ne2I1Os1Tl1() -> Molecule:
    """
    Load the molecule H6O2B2Ne2I1Os1Tl1 from 'fixtures/H6O2B2Ne2I1Os1Tl1.xyz'.
    """
    mol = Molecule("H6O2B2Ne2I1Os1Tl1")
    # get the Path of this file
    testsdir = Path(__file__).resolve().parents[1]
    xyz_file = testsdir / "fixtures/H6O2B2Ne2I1Os1Tl1.xyz"
    mol.read_xyz_from_file(xyz_file)
    return mol


@pytest.fixture
def mol_H2O2B2I1Os1() -> Molecule:
    """
    Load the molecule H2O2B2I1Os1 from 'fixtures/H2O2B2I1Os1.xyz'.
    """
    mol = Molecule("H2O2B2I1Os1")
    # get the Path of this file
    testsdir = Path(__file__).resolve().parents[1]
    xyz_file = testsdir / "fixtures/H2O2B2I1Os1.xyz"
    mol.read_xyz_from_file(xyz_file)
    return mol


@pytest.fixture
def mol_H2O2B2I1Os1_opt() -> Molecule:
    """
    Load the optimized molecule H2O2B2I1Os1 from 'fixtures/H2O2B2I1Os1_opt.xyz'.
    """
    mol = Molecule("H2O2B2I1Os1_opt")
    # get the Path of this file
    testsdir = Path(__file__).resolve().parents[1]
    xyz_file = testsdir / "fixtures/H2O2B2I1Os1_opt.xyz"
    mol.read_xyz_from_file(xyz_file)
    return mol


def test_detect_fragments_H6O2B2Ne2I1Os1Tl1(
    mol_H6O2B2Ne2I1Os1Tl1: Molecule,
    mol_H2O2B2I1Os1: Molecule,
) -> None:
    """
    Test the detection of fragments in the molecule H2O2B2I1Os1.
    """
    fragmols = detect_fragments(mol_H6O2B2Ne2I1Os1Tl1, verbosity=0)
    assert len(fragmols) == 7

    # check that the first fragment is equal to the molecule H2O2B2I1Os1
    assert np.allclose(fragmols[0].xyz, mol_H2O2B2I1Os1.xyz)


@pytest.mark.optional
def test_iterative_optimization(
    mol_H6O2B2Ne2I1Os1Tl1: Molecule, mol_H2O2B2I1Os1_opt: Molecule
) -> None:
    """
    Test the iterative optimization of the molecule H6O2B2Ne2I1Os1Tl1.
    """
    # initialize a configuration object
    config = ConfigManager()
    config.refine.max_frag_cycles = 1
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
    mol = mol_H6O2B2Ne2I1Os1Tl1
    mol.charge = 0
    mol_opt = iterative_optimization(
        mol,
        engine=engine,
        config_generate=config.generate,
        config_refine=config.refine,
        verbosity=0,
    )
    mol_ref = mol_H2O2B2I1Os1_opt

    # assert number of atoms in mol_opt is equal to number of atoms in mol_ref
    assert mol_opt.num_atoms == mol_ref.num_atoms
    # assert that the coordinates of mol_opt are close to the coordinates of mol_ref
    assert np.allclose(mol_opt.xyz, mol_ref.xyz, atol=1e-4)
