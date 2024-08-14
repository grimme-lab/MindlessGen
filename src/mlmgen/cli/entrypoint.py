"""
Entrypoint for command line interface.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import warnings

import toml

from ..generator import generator
from .cli_parser import cli_parser as cl


def console_entry_point(argv: Sequence[str] | None = None) -> int:
    """
    Entrypoint for command line interface.
    """
    # Step 1: Parse CLI arguments
    args = cl(argv)

    # Generate a default config that corresponds to this configuration file
    # after parsing the toml file into a dictionary
    DEFAULT_CONFIG = {
        "general": {
            "verbosity": 1,
            "engine": "xtb",
            "max_cycles": 100,
        },
        "xtb": {},
        "orca": {},
    }
    # Step 2: Find the configuration file (CLI provided or default search)
    config_file = find_config_file(args["general"]["config"])

    # Step 3: Load the configuration
    if config_file:
        print(f"Reading configuration from file: '{config_file}'")
        config = load_config(config_file)
    else:
        config = DEFAULT_CONFIG

    # Step 4: Merge with CLI arguments, giving precedence to CLI
    merged_config = merge_config_with_cli(config, args)

    # Use `final_config` in your program
    if merged_config["general"]["verbosity"] > 1:
        print(merged_config)
    raise SystemExit(generator(merged_config))


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


def load_config(config_file):
    """
    Load the configuration from the provided TOML file.
    """
    return toml.load(config_file)


def merge_config_with_cli(config, cli_args_dict):
    """
    Merge CLI arguments with the configuration, giving precedence to CLI.
    """

    for subcommand in cli_args_dict.keys():
        for key, value in cli_args_dict[subcommand].items():
            if value is not None:
                config[subcommand][key] = value
    return config
