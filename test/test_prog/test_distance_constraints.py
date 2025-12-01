from __future__ import annotations

import pytest

from mindlessgen.prog import DistanceConstraint, XTBConfig, ConfigManager  # type: ignore


def test_distance_constraint_cli_parsing():
    constraint = DistanceConstraint.from_cli_string("O,H,1.05")
    assert constraint.element_a == "O"
    assert constraint.element_b == "H"
    assert constraint.atomic_numbers == (8, 1)
    assert constraint.distance == pytest.approx(1.05)


def test_distance_constraint_from_mapping():
    data = {"pair": ["C", "N"], "distance": 2.0}
    constraint = DistanceConstraint.from_mapping(data)
    assert constraint.element_a == "C"
    assert constraint.element_b == "N"
    assert constraint.atomic_numbers == (6, 7)


def test_xtb_config_accepts_constraint_sequence():
    cfg = XTBConfig()
    cfg.distance_constraints = [{"pair": ["O", "H"], "distance": 1.0}]
    assert len(cfg.distance_constraints) == 1
    assert cfg.distance_constraints[0].element_a == "O"


def test_force_constant_validation():
    cfg = XTBConfig()
    cfg.distance_constraint_force_constant = 0.8
    assert cfg.distance_constraint_force_constant == pytest.approx(0.8)
    with pytest.raises(ValueError):
        cfg.distance_constraint_force_constant = 0


def test_check_config_allows_non_fixed_composition():
    config = ConfigManager()
    config.general.parallel = 1
    config.refine.ncores = 1
    config.xtb.distance_constraints = [DistanceConstraint.from_cli_string("He,He,2.0")]
    config.generate.fixed_composition = False
    config.generate.element_composition = "He:2-*"
    config.refine.engine = "xtb"

    config.check_config()


def test_check_config_requires_matching_counts_when_fixed():
    config = ConfigManager()
    config.general.parallel = 1
    config.refine.ncores = 1
    config.xtb.distance_constraints = [DistanceConstraint.from_cli_string("He,He,2.0")]
    config.generate.fixed_composition = True
    config.generate.element_composition = "He:1-1"
    config.refine.engine = "xtb"

    with pytest.raises(ValueError):
        config.check_config()


def test_check_config_accepts_matching_counts():
    config = ConfigManager()
    config.general.parallel = 1
    config.refine.ncores = 1
    config.xtb.distance_constraints = [DistanceConstraint.from_cli_string("O,H,1.0")]
    config.generate.fixed_composition = True
    config.generate.element_composition = "O:1-1, H:1-1"
    config.generate.min_num_atoms = 2
    config.generate.max_num_atoms = 2
    config.refine.engine = "xtb"

    config.check_config()
