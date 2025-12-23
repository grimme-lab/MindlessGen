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
            use_xtb_driver=False,
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


@pytest.fixture
def make_orca():
    def _factory(cfg=None, xtb_cfg_param=None):
        cfg = cfg or DummyORCAConfig()
        return ORCA(path="/usr/bin/orca", orcacfg=cfg, xtb_config=xtb_cfg_param)

    return _factory


def test_run_xtb_driver_success(monkeypatch, tmp_path, make_orca):
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
        str(Path("/fake/xtb")),
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


def test_run_xtb_driver_failure_returns_error(monkeypatch, tmp_path, make_orca):
    """Ensure the ORCA wrapper surfaces errors from the xTB driver."""
    orca = make_orca()
    monkeypatch.setattr(orca, "_get_xtb_executable", lambda: Path("/fake/xtb"))

    def fake_run(*_, **kwargs):
        del kwargs
        raise sp.CalledProcessError(1, "xtb", output=b"bad", stderr=b"worse")

    monkeypatch.setattr(sp, "run", fake_run)
    out, err, code = orca._run_xtb_driver(  # pylint: disable=protected-access
        tmp_path, "geom.xyz", "ctrl.inp", ncores=1
    )
    assert (out, err, code) == ("bad", "worse", 1)


def test_get_xtb_executable_raises_when_missing(monkeypatch, make_orca):
    orca = make_orca()

    def fake_get_xtb_path(candidate):
        raise ImportError("not found")

    monkeypatch.setattr("mindlessgen.qm.orca.get_xtb_path", fake_get_xtb_path)
    with pytest.raises(RuntimeError, match="xTB executable not found"):
        orca._get_xtb_executable()


def test_get_xtb_executable_prefers_xtb_cfg_path(monkeypatch, make_orca):
    xtb_cfg_constraints = DummyXTBConfig()
    xtb_cfg_constraints.xtb_path = "xtb_from_xtb_cfg"
    orca = make_orca(xtb_cfg_param=xtb_cfg_constraints)
    called = {}

    def fake_get_xtb_path(candidate):
        called.setdefault("candidates", []).append(candidate)
        return Path("/resolved/xtb_cfg")

    monkeypatch.setattr("mindlessgen.qm.orca.get_xtb_path", fake_get_xtb_path)
    assert orca._get_xtb_executable() == Path("/resolved/xtb_cfg")
    assert called["candidates"][0] == "xtb_from_xtb_cfg"


def test_should_use_xtb_driver_checks_distance_constraints(make_orca):
    cfg = DummyORCAConfig(use_xtb_driver=True)
    xtb_constraints = DummyXTBConfig(distance_constraints=[object()])
    orca = make_orca(cfg=cfg, xtb_cfg_param=xtb_constraints)
    assert orca._should_use_xtb_driver() is True

    no_constraints = DummyXTBConfig(distance_constraints=[])
    orca_no_constraints = make_orca(cfg=cfg, xtb_cfg_param=no_constraints)
    assert orca_no_constraints._should_use_xtb_driver() is False

    cfg_disabled = DummyORCAConfig(use_xtb_driver=False)
    orca_disabled = make_orca(cfg=cfg_disabled, xtb_cfg_param=xtb_constraints)
    assert orca_disabled._should_use_xtb_driver() is False
    orca_missing_xtb = make_orca(cfg=cfg, xtb_cfg_param=None)
    assert orca_missing_xtb._should_use_xtb_driver() is False


def test_write_xtb_input_creates_expected_file(monkeypatch, tmp_path, make_orca):
    xtb_cfg_instance = DummyXTBConfig(
        distance_constraints=["dummy"], distance_constraint_force_constant=0.7
    )
    orca = make_orca(xtb_cfg_param=xtb_cfg_instance)
    monkeypatch.setattr(orca, "_get_xtb_executable", lambda: Path("/fake/xtb"))

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
    orca._write_xtb_input(DummyMolecule(), target, "orca.inp")
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
