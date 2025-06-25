"""
Python script that is based on MindlessGen
and filters compounds based on allowed and forbidden elements.
"""

import argparse
from pathlib import Path
from tqdm import tqdm
from mindlessgen.molecules import Molecule  # type: ignore
from mindlessgen.molecules import get_molecules_from_filesystem  # type: ignore


def get_args() -> argparse.Namespace:
    """
    Get the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Detect fragments for a given list of molecules."
    )
    parser.add_argument(
        "--verbosity", "-v", type=int, default=1, help="Verbosity level."
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
        required=False,
        default=None,
        help="Allowed elements for the molecules. "
        + "Format example: `--allowed-elements '57-71, 81-*'",
    )
    parser.add_argument(
        "--required-elements-all",
        type=str,
        required=False,
        default=None,
        help="Required element(s) that MUST be in each molecule (ALL of them must be contained). "
        + "Format example: `--required-elements-all '57-71, 81-*'",
    )
    parser.add_argument(
        "--required-elements-one",
        type=str,
        required=False,
        default=None,
        help="Required element(s) that MUST be in each molecule "
        + "(at least one of them must be contained). "
        + "Format example: `--required-elements-one '57-71, 81-*'",
    )
    parser.add_argument(
        "--min-charge",
        type=int,
        required=False,
        default=None,
        help="Minimum charge for the molecules." + "Format example: `--min-charge -1`",
    )
    parser.add_argument(
        "--max-charge",
        type=int,
        required=False,
        default=None,
        help="Maximum charge for the molecules." + "Format example: `--max-charge 2`",
    )
    parser.add_argument(
        "--max-uhf",
        type=int,
        required=False,
        default=None,
        help="Maximum number of unpaired electrons (UHF) for the molecules."
        + " Format example: `--max-uhf 2`",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=False,
        default="selected_elements_molecules.list",
        help="Output file for the selected elements.",
    )
    return parser.parse_args()


def parse_element_list(allowed_elements: str) -> list[int]:
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


def molecule_has_required_elements(
    mol: Molecule, required_elements: list[tuple], verbosity: int
) -> bool:
    """
    Check whether a molecule contains the required elements.
    """
    # loop over all tuples of required element combinations
    contained_combinations: list[bool] = [False] * len(required_elements)
    for k, req_elem in enumerate(required_elements):
        # list of boolean values with the same length as the number of req_elem
        contained: list[bool] = [False] * len(req_elem)
        for i, ati in enumerate(req_elem):
            # check if the required element is in the molecule
            if ati in mol.ati:
                contained[i] = True
        # check if all elements of the respective required element combination are found
        if all(contained):
            contained_combinations[k] = True
    # check if any of the combinations is True
    if any(contained_combinations):
        if verbosity > 1:
            print(f"Molecule {mol.name} has the required elements.")
        return True
    if verbosity > 1:
        print(f"Molecule {mol.name} does not have the required elements.")
    return False


def main() -> int:
    """
    Main function that is called when the script is executed
    from the command line.
    """
    args = get_args()
    if (
        not args.allowed_elements
        and not args.required_elements_all
        and not args.required_elements_one
        and not args.min_charge
        and not args.max_charge
        and not args.max_uhf
    ):
        raise ValueError(
            "Either --allowed-elements, --required-elements_XXX, --min-charge, "
            + "--max-charge, or --max-uhf must be provided."
        )
    if args.required_elements_all and args.required_elements_one:
        raise ValueError(
            "Both --required-elements-all and "
            + "--required-elements-one cannot be provided at the same time."
        )
    if args.allowed_elements:
        allowed_elements = parse_element_list(args.allowed_elements)
    if args.required_elements_all:
        required_elements_all = parse_element_list(args.required_elements_all)
    if args.required_elements_one:
        required_elements_one = parse_element_list(args.required_elements_one)

    output_file = Path(args.output_file).resolve()
    if args.verbosity > 0:
        if args.allowed_elements:
            print(f"Allowed elements: {allowed_elements}")
        print(f"Output file: {output_file}")

    # required elements is a list of tuples
    # one tuple per set of required elements that must be contained at the same time
    # e.g. [(55, 56)] means that both 55 and 56 must be contained in the molecule
    # [(54),(55)] means that either 54 or 55 must be contained in the molecule
    required_elements: list[tuple] = []
    if args.required_elements_all:
        required_elements.append(tuple(required_elements_all))
    if args.required_elements_one:
        for elem in required_elements_one:
            required_elements.append(tuple([elem]))
    if args.verbosity > 0:
        print(f"Required elements: {required_elements}")

    mols = get_molecules_from_filesystem(keyword=args.keyword, verbosity=args.verbosity)
    with open(output_file, "w", encoding="utf8") as sel_elem_file:
        for mol in tqdm(mols, desc="Checking composition...", unit="molecule"):
            # check if all elements in the molecule are allowed
            if args.allowed_elements:
                if all(ati in allowed_elements for ati in mol.ati):
                    if args.verbosity > 1:
                        print(f"Molecule {mol.name} has only allowed elements.")
                else:
                    if args.verbosity > 1:
                        print(f"Molecule {mol.name} has forbidden elements.")
                    continue
            if required_elements and (
                not molecule_has_required_elements(
                    mol, required_elements, args.verbosity
                )
            ):
                continue

            if args.min_charge is not None and mol.charge < args.min_charge:
                if args.verbosity > 1:
                    print(f"Molecule {mol.name} has charge {mol.charge}.")
                continue
            if args.max_charge is not None and mol.charge > args.max_charge:
                if args.verbosity > 1:
                    print(f"Molecule {mol.name} has charge {mol.charge}.")
                continue
            if args.max_uhf is not None and mol.uhf > args.max_uhf:
                if args.verbosity > 1:
                    print(f"Molecule {mol.name} has UHF {mol.uhf}.")
                continue

            sel_elem_file.write(mol.name + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
