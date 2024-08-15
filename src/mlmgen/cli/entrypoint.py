"""
Entrypoint for command line interface.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import warnings


from ..generator import generator
from .cli_parser import cli_parser as cl
from ..prog import ConfigManager
from ..molecules import Molecule


def console_entry_point(argv: Sequence[str] | None = None) -> int:
    """
    Entrypoint for command line interface.
    """
    # Step 1: Parse CLI arguments
    args = cl(argv)

    # Step 2: Find the configuration file (CLI provided or default search)
    config_file = find_config_file(args["general"]["config"])

    # Step 3: Load the configuration
    if config_file:
        print(f"Reading configuration from file: '{config_file}'")
        config = ConfigManager(config_file)
    else:
        config = ConfigManager()

    # Step 4: Merge with CLI arguments, giving precedence to CLI
    config.load_from_dict(args)

    # Step 5: Run the generator
    if config.general.verbosity > 1:
        print(config)

    try:
        returnobject: Molecule | int | None = generator(config)
    except RuntimeError as e:
        print(f"\nGeneration failed: {e}")
        raise RuntimeError("Generation failed.") from e
    if isinstance(returnobject, int):
        raise SystemExit(returnobject)
    if isinstance(returnobject, Molecule):
        returnobject.write_xyz_to_file()
        raise SystemExit(0)
    raise SystemExit(1)


def find_config_file(cli_config_path: str | Path | None = None) -> Path | None:
    """
    Finds the configuration file. If a path is provided via CLI, use it.
    Otherwise, search in predefined locations.
    """
    # CLI provided config file
    if cli_config_path:
        config_path = Path(cli_config_path).resolve()
        if config_path.is_file():
            return config_path
        raise FileNotFoundError(f"Configuration file not found at {cli_config_path}")

    # Search paths
    search_paths = [
        Path.home() / "mlmgen.toml",  # $USER/mlmgen.toml
        Path.cwd() / "mlmgen.toml",  # Current directory
    ]

    # Find the config file
    for path in search_paths:
        if path.is_file():
            return path

    # If no config file is found, raise a warning
    warnings.warn("No configuration file found. Using default configuration.")
    return None
