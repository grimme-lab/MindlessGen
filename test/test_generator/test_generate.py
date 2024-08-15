import pytest

from mlmgen.generator import generator  # type: ignore
from mlmgen.prog import ConfigManager  # type: ignore
from mlmgen.molecules import Molecule  # type: ignore


@pytest.mark.optional
def test_generator():
    config = ConfigManager()
    config.general.engine = "xtb"
    config.general.max_cycles = 10000
    config.general.parallel = 8
    config.general.min_num_atoms = 2
    config.general.max_num_atoms = 100
    config.general.verbosity = 0

    molecule, exitcode = generator(config)
    assert exitcode == 0
    assert isinstance(molecule, Molecule)
