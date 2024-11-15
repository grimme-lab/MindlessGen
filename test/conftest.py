"""
Setup for pytest.
"""

import numpy as np
import pytest

np.set_printoptions(precision=16)


def pytest_addoption(parser):
    parser.addoption(
        "--optional",
        action="store_true",
        default=False,
        help="Run optional tests",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--optional"):
        # If --optional is provided, do nothing, so all tests run
        return
    skip_optional = pytest.mark.skip(reason="Need '--optional' option to run")
    for item in items:
        if "optional" in item.keywords:
            item.add_marker(skip_optional)
