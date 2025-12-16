import subprocess as sp
from pathlib import Path
from types import SimpleNamespace
import pytest
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
        )
        defaults.update(kwargs)
        super().__init__(**defaults)


class DummyXTBConfig(SimpleNamespace):
    def __init__(self, **kwargs):
        defaults = dict(
            distance_constraints=None, distance_constraint_force_constant=None
        )
        defaults.update(kwargs)
        super().__init__(**defaults)


class DummyMolecule:
    ati = [1, 1, 6, 8]


def make_orca(cfg=None, xtb_cfg=None):
    cfg = cfg or DummyORCAConfig()
    return ORCA(path="/usr/bin/orca", orcacfg=cfg, xtb_config=xtb_cfg)


def test_run_xtb_driver_success(monkeypatch, tmp_path):
    orca = make_orca(cfg=DummyORCAConfig(optlevel="tight"))
    monkeypatch.setattr(orca, "_get_xtb_executable", lambda: Path("/fake/xtb"))
    captured = {}

    def fake_run(args, cwd, capture_output, check):
        captured["args"] = args
        assert cwd == tmp_path
        assert capture_output and check
        return SimpleNamespace(stdout=b"ok", stderr=b"")

    monkeypatch.setattr(sp, "run", fake_run)
    out, err, code = orca._run_xtb_driver(tmp_path, "geom.xyz", "ctrl.inp", ncores=4)
    assert captured["args"] == [
        "/fake/xtb",
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


def test_run_xtb_driver_failure_returns_error(monkeypatch, tmp_path):
    orca = make_orca()
    monkeypatch.setattr(orca, "_get_xtb_executable", lambda: Path("/fake/xtb"))

    def fake_run(*_, **__):
        raise sp.CalledProcessError(1, "xtb", output=b"bad", stderr=b"worse")

    monkeypatch.setattr(sp, "run", fake_run)
    out, err, code = orca._run_xtb_driver(  # pylint: disable=protected-access
        tmp_path, "geom.xyz", "ctrl.inp", ncores=1
    )
    assert (out, err, code) == ("bad", "worse", 1)


def test_get_xtb_executable_prefers_configured_path(monkeypatch):
    cfg = DummyORCAConfig(xtb_driver_path="custom_xtb")
    orca = make_orca(cfg=cfg)
    called = {}

    def fake_get_xtb_path(candidate):
        called["candidate"] = candidate
        return Path("/resolved/xtb")

    monkeypatch.setattr("mindlessgen.qm.orca.get_xtb_path", fake_get_xtb_path)
    assert orca._get_xtb_executable() == Path("/resolved/xtb")
    assert called["candidate"] == "custom_xtb"


def test_get_xtb_executable_raises_when_missing(monkeypatch):
    orca = make_orca()

    def fake_get_xtb_path(candidate):
        raise ImportError("not found")

    monkeypatch.setattr("mindlessgen.qm.orca.get_xtb_path", fake_get_xtb_path)
    with pytest.raises(RuntimeError, match="xTB executable not found"):
        orca._get_xtb_executable()


def test_should_use_xtb_driver_checks_distance_constraints():
    orca = make_orca(xtb_cfg=DummyXTBConfig(distance_constraints=[object()]))
    assert orca._should_use_xtb_driver() is True
    orca_no_constraints = make_orca(xtb_cfg=DummyXTBConfig(distance_constraints=[]))
    assert orca_no_constraints._should_use_xtb_driver() is False


def test_write_xtb_input_creates_expected_file(monkeypatch, tmp_path):
    xtb_cfg = DummyXTBConfig(
        distance_constraints=["dummy"], distance_constraint_force_constant=0.7
    )
    orca = make_orca(xtb_cfg=xtb_cfg)
    monkeypatch.setattr(
        ORCA,
        "_prepare_distance_constraint_section",
        lambda self, mol: ["  distance: 1, 2, 1.00000"],
    )
    target = tmp_path / "xtb.inp"
    orca._write_xtb_input(DummyMolecule(), target, "orca.inp")
    content = target.read_text().splitlines()
    assert content[:4] == [
        "$constrain",
        "  force constant= 0.7",
        "  distance: 1, 2, 1.00000",
        "$end",
    ]
    assert "$external" in content
    assert "  orca input file= orca.inp" in content
    assert f"  orca bin= {orca.path}" in content
