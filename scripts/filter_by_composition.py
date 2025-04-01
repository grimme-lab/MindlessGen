"""
Python script that is based on MindlessGen
and filters compounds based on allowed and forbidden elements.
"""

import argparse
from pathlib import Path
from tqdm import tqdm
from mindlessgen.molecules import Molecule  # type: ignore


def get_molecules_from_filesystem(keyword: str) -> list[Molecule]:
    """
    Get a list of molecules from the filesystem.
    """
    # check if the file exists
    if not Path(keyword).exists():
        raise FileNotFoundError(f"File '{keyword}' does not exist.")
    # read the file
    with open(keyword, encoding="utf-8") as file:
        mol_names = file.readlines()
    # get the molecules and return them
    mol_list: list[Molecule] = []
    for mol_name in tqdm(
        mol_names, desc="Processing molecules from files...", unit="molecule"
    ):
        mol_name = mol_name.strip()
        mol = Molecule.read_mol_from_file(mol_name + ".xyz")
        mol_list.append(mol)
    return mol_list


def get_args() -> argparse.Namespace:
    """
    Get the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Detect fragments for a given list of molecules."
    )
    parser.add_argument(
        "--keyword",
        type=str,
        required=False,
        default="molecules.list",
        help="Keyword for the file that contains the list of molecules.",
    )
    parser.add_argument(
        "--allowed-elements",
        type=str,
        required=True,
        default=None,
        help="Allowed elements for the molecules. Format example: `--allowed-elements '57-71, 81-*'",
    )
    return parser.parse_args()


def parse_allowed_elements(allowed_elements: str) -> list[int]:
    """
    Parse the allowed elements from a string.
    """
    set_allowed_elements: set[int] = set()
    elements = allowed_elements.split(",")
    elements = [element.strip() for element in elements]

    for item in elements:
        if "-" in item:
            start: str | int
            end: str | int
            start, end = item.split("-")
            if end == "*" and start == "*":
                raise ValueError("Both start and end cannot be wildcard '*'.")
            if end == "*":
                end = 103  # Set to a the maximum atomic number in mindlessgen
            if start == "*":
                start = 0
            set_allowed_elements.update(
                range(int(start) - 1, int(end))
            )  # Subtract 1 to convert to 0-based indexing
        else:
            set_allowed_elements.add(
                int(item) - 1
            )  # Subtract 1 to convert to 0-based indexing

    return sorted(list(set_allowed_elements))


def main() -> int:
    """
    Main function that is called when the script is executed
    from the command line.
    """
    args = get_args()
    mols = get_molecules_from_filesystem(keyword=args.keyword)
    allowed_elements = parse_allowed_elements(args.allowed_elements)
    with open(
        "selected_elements_molecules.list", "w", encoding="utf8"
    ) as sel_elem_file:
        for mol in tqdm(mols, desc="Detecting fragments...", unit="molecule"):
            # check if all elements in the molecule are allowed
            if all(ati in allowed_elements for ati in mol.ati):
                print(f"Molecule {mol.name} has only allowed elements.")
                sel_elem_file.write(mol.name + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
