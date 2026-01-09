"""
Python script that is based on MindlessGen
and filters compounds that are redundant stereoisomers.
"""

import argparse
from pathlib import Path
from collections import defaultdict

from tqdm import tqdm
import networkx as nx  # type: ignore

from mindlessgen.molecules import get_molecules_from_filesystem  # type: ignore
from mindlessgen.molecules import get_molecular_graph  # type: ignore


def get_args() -> argparse.Namespace:
    """
    Get the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Detect stereoisomers for a given list of molecules."
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
        "--output-file",
        type=str,
        required=False,
        default="selected_elements_molecules.list",
        help="Output file for the selected elements.",
    )
    return parser.parse_args()


def main() -> int:
    """
    Main function that is called when the script is executed.
    """
    args = get_args()
    output_file = Path(args.output_file).resolve()
    if args.verbosity > 0:
        print(f"Output file: {output_file}")
    mols = get_molecules_from_filesystem(keyword=args.keyword, verbosity=args.verbosity)

    seen_hashes: defaultdict[str, list[str]] = defaultdict(
        list
    )  # maps graph hashes to list of mol indices or names

    with open(output_file, "w", encoding="utf8") as sel_elem_file:
        for i, mol in enumerate(
            tqdm(mols, desc="Checking composition...", unit="molecule")
        ):
            graph = get_molecular_graph(mol, 1.25, verbosity=args.verbosity)

            # Get WL hash with atom type info
            g_hash = nx.weisfeiler_lehman_graph_hash(graph, node_attr="element")

            if g_hash in seen_hashes.keys():
                if args.verbosity > 1:
                    print(
                        f"Found stereoisomer: {seen_hashes[g_hash]} "
                        + f"and {mol.name} with hash {g_hash}"
                    )
                seen_hashes[g_hash].append(mol.name)
                continue
            seen_hashes[g_hash].append(mol.name)
            sel_elem_file.write(mol.name + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
