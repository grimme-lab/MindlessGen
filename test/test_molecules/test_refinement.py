"""
Test the postprocessing functions in the molecules module.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest
from mindlessgen.prog import ConfigManager, GenerateConfig  # type: ignore
from mindlessgen.molecules import detect_fragments  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore
from mindlessgen.molecules import iterative_optimization  # type: ignore
from mindlessgen.qm import XTB, get_xtb_path  # type: ignore

TESTSDIR = Path(__file__).resolve().parents[1]


@pytest.fixture
def mol_C13H14() -> Molecule:
    """
    Load the molecule C13H14 from 'fixtures/C13H14.xyz'.
    """
    molfile = TESTSDIR / "fixtures/C13H14.xyz"
    mol = Molecule.read_mol_from_file(molfile)
    return mol


@pytest.fixture
def mol_C7H8() -> Molecule:
    """
    Load the molecule C7H8 from 'fixtures/C7H8.xyz'.
    """
    molfile = TESTSDIR / "fixtures/C7H8.xyz"
    mol = Molecule.read_mol_from_file(molfile)
    return mol


@pytest.fixture
def mol_C2H2N2O1Co1Ge2Ta1Hg1() -> Molecule:
    """
    Load the molecule C2H2N2O1Co1Ge2Ta1Hg1 from 'fixtures/C2H2N2O1Co1Ge2Ta1Hg1.xyz'.
    """
    testsdir = Path(__file__).resolve().parents[1]
    molfile = testsdir / "fixtures/C2H2N2O1Co1Ge2Ta1Hg1.xyz"
    mol = Molecule.read_mol_from_file(molfile)
    return mol


@pytest.fixture
def mol_C2H1N2O1Co1Ge2Ta1Hg1_frag() -> Molecule:
    """
    Load the optimized molecule C2H1N2O1Co1Ge2Ta1Hg1 from 'fixtures/C2H1N2O1Co1Ge2Ta1Hg1_frag.xyz'.
    """
    testsdir = Path(__file__).resolve().parents[1]
    molfile = testsdir / "fixtures/C2H1N2O1Co1Ge2Ta1Hg1_frag.xyz"
    mol = Molecule.read_mol_from_file(molfile)
    return mol


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
    mol.charge = 0
    mol.uhf = 0
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
    fragmols = detect_fragments(
        mol=mol_H6O2B2Ne2I1Os1Tl1,
        molecular_charge=GenerateConfig().molecular_charge,
        vdw_scaling=1.3333,
        verbosity=0,
    )
    assert len(fragmols) == 7

    # check that the first fragment is equal to the molecule H2O2B2I1Os1
    assert np.allclose(fragmols[0].xyz, mol_H2O2B2I1Os1.xyz)


@pytest.mark.optional
def test_iterative_optimization(mol_C13H14: Molecule, mol_C7H8: Molecule) -> None:
    """
    Test the iterative optimization of the molecule H6O2B2Ne2I1Os1Tl1.
    """
    # initialize a configuration object
    config = ConfigManager()
    config.generate.min_num_atoms = 2
    config.generate.max_num_atoms = 100
    config.refine.hlgap = 0.001  # TODO: Change charge assignment such that
    # fragment charge is not completely random anymore. Currently, that's the
    # reason for a virtually switched off HL gap check (fragment can be -2, 0, 2)
    config.refine.max_frag_cycles = 1
    if config.refine.engine == "xtb":
        try:
            xtb_path = get_xtb_path()
            if not xtb_path:
                raise ImportError("xtb not found.")
        except ImportError as e:
            raise ImportError("xtb not found.") from e
        engine = XTB(xtb_path, config.xtb)
    else:
        raise NotImplementedError("Engine not implemented.")
    mol = mol_C13H14
    mol_opt = iterative_optimization(
        mol,
        engine=engine,
        config_generate=config.generate,
        config_refine=config.refine,
        verbosity=2,
    )
    mol_ref = mol_C7H8

    # assert number of atoms in mol_opt is equal to number of atoms in mol_ref
    assert mol_opt.num_atoms == mol_ref.num_atoms
    # assert that the coordinates of mol_opt are close to the coordinates of mol_ref
    assert np.allclose(mol_opt.xyz, mol_ref.xyz, atol=1e-4)
