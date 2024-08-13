from pathlib import Path
import networkx as nx
import numpy as np
from ..qm.base import QMMethod
from .molecule import Molecule
from .miscellaneous import set_random_charge


def postprocess(mol: Molecule, engine: QMMethod, verbosity: int = 1) -> Molecule:
    """
    Postprocess the molecule.
    """
    optmol = engine.optimize(mol)
    fragmols = detect_fragments(optmol, verbosity)

    if verbosity > 1:
        for i, fragmol in enumerate(fragmols):
            # generate a directory for each fragment
            Path(f"fragment_{i}").mkdir(exist_ok=True)
            fragmol.write_xyz_to_file(f"fragment_{i}/fragment_{i}.xyz")

    # return first fragment
    if verbosity > 1:
        print(f"Returning fragment 0: {fragmols[0]}")
    return fragmols[0]


def detect_fragments(mol: Molecule, verbosity: int = 1) -> list[Molecule]:
    """
    Detects fragments in a molecular system based on atom positions and covalent radii.

    CAUTION: 0-based indexing is used for atoms and molecules!

    Parameters:
    - mol: A Molecule object, containg the xyz coordinates and the types of atoms.
    - covalent_radii: A list of covalent radii corresponding to each atom.

    Returns:
    - A list of fragments, where each fragment is a Molecule object.
    """
    graph = nx.Graph()

    # Add nodes (atoms) to the graph
    graph.add_nodes_from(range(mol.num_atoms))

    # Calculate pairwise distances and add edges if atoms are bonded
    for i in range(mol.num_atoms - 1):
        for j in range(i + 1, mol.num_atoms):
            distance = np.linalg.norm(mol.xyz[i] - mol.xyz[j])
            sum_radii = get_cov_radii(mol.ati[i]) + get_cov_radii(mol.ati[j])
            if distance <= sum_radii * 0.9:
                graph.add_edge(i, j)

    # Detect connected components (fragments)
    fragments = [list(component) for component in nx.connected_components(graph)]
    # sort the fragments by len
    fragments = sorted(fragments, key=len, reverse=True)

    # Generate a Molecule object for each fragment
    fragment_molecules = []
    for fragment in fragments:
        if verbosity > 1:
            print(f"Fragment: {fragment}")
        fragment_molecule = mol.copy()
        # Update the number of atoms
        fragment_molecule.num_atoms = len(fragment)
        # Remove the atoms that are not in the fragment from xyz
        fragment_molecule.xyz = np.array([mol.xyz[i] for i in fragment])
        # Remove the atoms that are not in the fragment from ati
        fragment_molecule.ati = np.array([mol.ati[i] for i in fragment])
        # Update atlist by going through ati and add the occurences of each atom onto an empty array
        fragment_molecule.atlist = np.zeros(102, dtype=int)
        for atom in fragment_molecule.ati:
            fragment_molecule.atlist[atom] += 1
        # Update the charge of the fragment molecule
        fragment_molecule.charge = set_random_charge(fragment_molecule.ati)
        if verbosity > 1:
            print(f"Fragment molecule: {fragment_molecule}")
        # Append the fragment molecule to the list
        fragment_molecules.append(fragment_molecule)

    return fragment_molecules


def get_cov_radii(at: int) -> float:
    """
    D3 covalent radii.
    """
    rcov = [
        0.80628308,
        1.15903197,
        3.02356173,
        2.36845659,
        1.94011865,
        1.88972601,
        1.78894056,
        1.58736983,
        1.61256616,
        1.68815527,
        3.52748848,
        3.14954334,
        2.84718717,
        2.62041997,
        2.77159820,
        2.57002732,
        2.49443835,
        2.41884923,
        4.43455700,
        3.88023730,
        3.35111422,
        3.07395437,
        3.04875805,
        2.77159820,
        2.69600923,
        2.62041997,
        2.51963467,
        2.49443835,
        2.54483100,
        2.74640188,
        2.82199085,
        2.74640188,
        2.89757982,
        2.77159820,
        2.87238349,
        2.94797246,
        4.76210950,
        4.20778980,
        3.70386304,
        3.50229216,
        3.32591790,
        3.12434702,
        2.89757982,
        2.84718717,
        2.84718717,
        2.72120556,
        2.89757982,
        3.09915070,
        3.22513231,
        3.17473967,
        3.17473967,
        3.09915070,
        3.32591790,
        3.30072128,
        5.26603625,
        4.43455700,
        4.08180818,
        3.70386304,
        3.98102289,
        3.95582657,
        3.93062995,
        3.90543362,
        3.80464833,
        3.82984466,
        3.80464833,
        3.77945201,
        3.75425569,
        3.75425569,
        3.72905937,
        3.85504098,
        3.67866672,
        3.45189952,
        3.30072128,
        3.09915070,
        2.97316878,
        2.92277614,
        2.79679452,
        2.82199085,
        2.84718717,
        3.32591790,
        3.27552496,
        3.27552496,
        3.42670319,
        3.30072128,
        3.47709584,
        3.57788113,
        5.06446567,
        4.56053862,
        4.20778980,
        3.98102289,
        3.82984466,
        3.85504098,
        3.88023730,
        3.90543362,
    ]
    return rcov[at]
