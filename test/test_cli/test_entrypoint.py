"""
Test program call from command line.
"""

import pytest

from mlmgen.cli import console_entry_point


def test_entrypoint(capsys: pytest.CaptureFixture) -> None:
    # pylint: disable=too-many-function-args
    console_entry_point(["2.0"])

    out, err = capsys.readouterr()
    assert out == "4.0\n"
    assert err == ""
