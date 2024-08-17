import pytest

from mindlessgen.generator import generator  # type: ignore
from mindlessgen.prog import ConfigManager  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore


@pytest.mark.optional
def test_generator():
    config = ConfigManager()
    config.general.engine = "xtb"
    config.general.max_cycles = 10000
    config.general.parallel = 8
    config.general.verbosity = 0

    molecule, exitcode = generator(config)
    assert exitcode == 0
    assert isinstance(molecule, Molecule)
