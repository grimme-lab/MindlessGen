"""
Python script that is based on MindlessGen
and is used to detect fragments for a given list of molecules.
"""

import argparse
from pathlib import Path
from tqdm import tqdm
from mindlessgen.molecules import Molecule, detect_fragments  # type: ignore


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
        "--min-num-atoms",
        type=int,
        required=False,
        default=7,
        help="Minimum number of atoms for a fragment to keep "
        + "(will be protocolled in a separate folder).",
    )
    return parser.parse_args()


def main() -> int:
    """
    Main function that is called when the script is executed
    from the command line.
    """
    args = get_args()
    mols = get_molecules_from_filesystem(keyword=args.keyword)
    # create new directory "new_single_molecules" if it does not exist
    newmoldir = Path("fragments").resolve()
    newmoldir.mkdir(exist_ok=True, parents=True)
    with open("nci_molecules.list", "w", encoding="utf8") as nci_file:
        with open("single_molecules.list", "w", encoding="utf8") as single_file:
            with open(
                newmoldir / "single_molecules.list", "w", encoding="utf8"
            ) as single_file_new_moldir:
                for mol in tqdm(mols, desc="Detecting fragments...", unit="molecule"):
                    fragmented_molecules = detect_fragments(
                        mol=mol,
                        molecular_charge=mol.charge,
                        vdw_scaling=1.25,
                        verbosity=1,
                    )
                    if len(fragmented_molecules) > 1:
                        nci_file.write(mol.name + "\n")
                        if fragmented_molecules[0].num_atoms >= args.min_num_atoms:
                            fragmented_molecules[0].write_xyz_to_file(
                                filename=(
                                    newmoldir
                                    / Path(
                                        "mlm_" + fragmented_molecules[0].name + ".xyz"
                                    )
                                ).resolve()
                            )
                            single_file_new_moldir.write(mol.name + "\n")
                    else:
                        single_file.write(mol.name + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
