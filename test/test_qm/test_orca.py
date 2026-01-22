import subprocess as sp
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from mindlessgen.molecules import Molecule  # type: ignore
from mindlessgen.prog import DistanceConstraint, XTBConfig  # type: ignore
from mindlessgen.qm.orca import ORCA


class DummyORCAConfig(SimpleNamespace):
    def __init__(self, **kwargs):
        defaults = dict(
            functional="B3LYP",
            basis="def2-SVP",
            gridsize=2,
            scf_cycles=50,
            optlevel="",
            xtb_driver_path=None,
            xtb_path=None,
            use_xtb_driver=False,
        )
        defaults.update(kwargs)
        super().__init__(**defaults)


@pytest.fixture
def make_orca():
    def _factory(cfg=None, xtb_cfg_param=None):
        cfg = cfg or DummyORCAConfig()
        return ORCA(path="/usr/bin/orca", orcacfg=cfg, xtb_config=xtb_cfg_param)

    return _factory


@pytest.fixture
def xtb_cfg_oh():
    """
    xTB config with an O-H distance constraint.
    """
    cfg = XTBConfig()
    cfg.distance_constraints = [DistanceConstraint.from_cli_string("O,H,1.0")]
    cfg.distance_constraint_force_constant = 0.5
    return cfg


@pytest.fixture
def xtb_cfg_ff():
    """
    xTB config with an F-F distance constraint.
    """
    cfg = XTBConfig()
    cfg.distance_constraints = [DistanceConstraint.from_cli_string("F,F,1.0")]
    return cfg


@pytest.fixture
def xtb_cfg_fe3():
    """
    xTB config with a Fe-Fe distance constraint.
    """
    cfg = XTBConfig()
    cfg.distance_constraints = [
        DistanceConstraint.from_mapping({"pair": ["Fe", "Fe"], "distance": 2.5})
    ]
    return cfg


@pytest.fixture
def mol_oh():
    """
    Simple O-H molecule for constraint tests.
    """
    mol = Molecule("OH")
    mol.ati = np.array([7, 0])
    return mol


@pytest.fixture
def mol_h2():
    """
    Simple H2 molecule for constraint tests.
    """
    mol = Molecule("H2")
    mol.ati = np.array([0, 0])
    return mol


@pytest.fixture
def mol_fe3():
    """
    Simple Fe3 molecule for constraint tests.
    """
    mol = Molecule("Fe3")
    mol.ati = np.array([25, 25, 25])
    return mol


@pytest.fixture
def fake_xtb_path(monkeypatch):
    """
    Force xTB path discovery to a fake binary.
    """
    xtb_path = Path("/fake/xtb")
    monkeypatch.setattr("mindlessgen.qm.orca.get_xtb_path", lambda: xtb_path)
    return xtb_path


def test_run_xtb_driver_success(monkeypatch, tmp_path, make_orca, fake_xtb_path):
    orca = make_orca(
        cfg=DummyORCAConfig(optlevel="tight"),
        xtb_cfg_param=XTBConfig(),
    )
    captured = {}

    def fake_run(args, cwd, capture_output, check):
        captured["args"] = args
        assert cwd == tmp_path
        assert capture_output and check
        return SimpleNamespace(stdout=b"ok", stderr=b"")

    monkeypatch.setattr("mindlessgen.qm.xtb.sp.run", fake_run)
    out, err, code = orca._run_xtb_driver(tmp_path, "geom.xyz", "ctrl.inp", ncores=4)
    assert captured["args"] == [
        str(fake_xtb_path),
        "geom.xyz",
        "--opt",
        "tight",
        "--orca",
        "-I",
        "ctrl.inp",
    ]
    assert out == "ok"
    assert err == ""
    assert code == 0


def test_run_xtb_driver_failure_returns_error(
    monkeypatch, tmp_path, make_orca, fake_xtb_path
):
    """Ensure the ORCA wrapper surfaces errors from the xTB driver."""
    orca = make_orca(xtb_cfg_param=XTBConfig())

    def fake_run(*_, **kwargs):
        del kwargs
        raise sp.CalledProcessError(1, "xtb", output=b"bad", stderr=b"worse")

    monkeypatch.setattr("mindlessgen.qm.xtb.sp.run", fake_run)
    out, err, code = orca._run_xtb_driver(  # pylint: disable=protected-access
        tmp_path, "geom.xyz", "ctrl.inp", ncores=1
    )
    assert (out, err, code) == ("bad", "worse", 1)


def test_run_xtb_driver_requires_xtb_cfg(
    monkeypatch, tmp_path, make_orca, fake_xtb_path
):
    orca = make_orca()
    with pytest.raises(RuntimeError, match="xTB driver requested"):
        orca._run_xtb_driver(tmp_path, "geom.xyz", "ctrl.inp", ncores=1)


def test_write_xtb_input_creates_expected_file(
    monkeypatch, tmp_path, make_orca, fake_xtb_path, mol_oh
):
    xtb_cfg_instance = XTBConfig()
    xtb_cfg_instance.distance_constraints = []
    xtb_cfg_instance.distance_constraint_force_constant = 0.7
    orca = make_orca(xtb_cfg_param=xtb_cfg_instance)

    def fake_prepare(self, molecule, temp_dir):
        assert temp_dir == tmp_path
        (temp_dir / "xtb.inp").write_text(
            "\n".join(
                [
                    "$constrain",
                    " force constant= 0.7",
                    " distance: 1, 2, 1.00000",
                    "$end",
                    "",
                ]
            ),
            encoding="utf8",
        )
        return True

    monkeypatch.setattr(
        "mindlessgen.qm.orca.XTB._prepare_distance_constraint_file", fake_prepare
    )
    target = tmp_path / "xtb.inp"
    orca._write_xtb_input(mol_oh, target, "orca.inp")
    content = target.read_text().splitlines()
    assert content[:4] == [
        "$constrain",
        " force constant= 0.7",
        " distance: 1, 2, 1.00000",
        "$end",
    ]
    assert "$external" in content
    assert "  orca input file= orca.inp" in content
    assert f"  orca bin= {orca.path}" in content


def test_write_xtb_input_generates_constraints(
    fake_xtb_path, tmp_path, make_orca, xtb_cfg_oh, mol_oh
):
    orca = make_orca(xtb_cfg_param=xtb_cfg_oh)

    xtb_input = tmp_path / "xtb.inp"
    orca._write_xtb_input(mol_oh, xtb_input, "orca.inp")

    contents = xtb_input.read_text(encoding="utf8").splitlines()
    assert contents[0] == "$constrain"
    assert "force constant= 0.5" in contents[1]
    distance_lines = [line for line in contents if line.startswith(" distance:")]
    assert distance_lines == [" distance: 1, 2, 1.0"]
    assert "$external" in contents
    assert "  orca input file= orca.inp" in contents
    assert f"  orca bin= {orca.path}" in contents


def test_write_xtb_input_missing_atoms(
    fake_xtb_path, tmp_path, make_orca, xtb_cfg_ff, mol_h2
):
    orca = make_orca(xtb_cfg_param=xtb_cfg_ff)

    with pytest.raises(RuntimeError):
        orca._write_xtb_input(mol_h2, tmp_path / "xtb.inp", "orca.inp")


def test_distance_constraints_use_first_atoms(
    fake_xtb_path, tmp_path, make_orca, xtb_cfg_fe3, mol_fe3
):
    orca = make_orca(xtb_cfg_param=xtb_cfg_fe3)

    xtb_input = tmp_path / "xtb.inp"
    orca._write_xtb_input(mol_fe3, xtb_input, "orca.inp")

    contents = xtb_input.read_text(encoding="utf8").splitlines()
    distance_lines = [line for line in contents if "distance:" in line]

    assert len(distance_lines) == 1
    assert "1, 2" in distance_lines[0]
    assert ", 3," not in distance_lines[0]
