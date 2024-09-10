import pytest

from mindlessgen.generator import generator  # type: ignore
from mindlessgen.prog import ConfigManager  # type: ignore
from mindlessgen.molecules import Molecule  # type: ignore


@pytest.mark.optional
def test_generator():
    config = ConfigManager()
    config.refine.engine = "xtb"
    config.general.max_cycles = 10000
    config.general.parallel = 4
    config.general.verbosity = -1
    config.general.postprocess = False
    config.general.write_xyz = False

    molecules, exitcode = generator(config)
    assert exitcode == 0
    for molecule in molecules:
        assert isinstance(molecule, Molecule)
