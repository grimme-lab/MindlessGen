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
        config = ConfigManager(config_file)
    else:
        config = ConfigManager()

    # Step 4: Merge with CLI arguments, giving precedence to CLI
    config.load_from_dict(args)
    if config.general.verbosity >= 0:
        print(f"Reading configuration from file: '{config_file}'")

    try:
        molecules, exitcode = generator(config)
    except RuntimeError as e:
        print(f"\nGeneration failed: {e}")
        raise RuntimeError("Generation failed.") from e

    if molecules:
        # <- (Some) molecules were generated: Entry point for any post-processing of the generated molecules
        if exitcode != 0:
            warnings.warn(
                "Generation completed with errors for parts of the generation."
            )
            raise SystemExit(exitcode)
    if exitcode == 0:
        raise SystemExit(0)
    warnings.warn("CAUTION: Generation failed for all molecules. No files written.")
    raise SystemExit(exitcode)


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
        Path.cwd() / "mindlessgen.toml",  # Current directory
        Path.home() / "mindlessgen.toml",  # $USER/mindlessgen.toml
    ]

    # Find the config file
    for path in search_paths:
        if path.is_file():
            return path

    # If no config file is found, raise a warning
    warnings.warn("No configuration file found. Using default configuration.")
    return None
