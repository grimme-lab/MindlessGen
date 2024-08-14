from pathlib import Path
import networkx as nx  # type: ignore
import numpy as np
from ..qm.base import QMMethod
from .molecule import Molecule
from .miscellaneous import set_random_charge

COV_RADII = "pyykko"


def postprocess(
    mol: Molecule,
    engine: QMMethod,
    min_nat: int = 1,
    max_nat: int = 100,
    verbosity: int = 1,
) -> Molecule:
    """
    Postprocess the molecule.
    """
    # Optimize the initial random molecule
    try:
        optmol = engine.optimize(mol, verbosity)
    except RuntimeError as e:
        raise RuntimeError(f"First optimization failed: {e}") from e
    # Get all fragments
    fragmols = detect_fragments(optmol, verbosity)

    if verbosity > 1:
        for i, fragmol in enumerate(fragmols):
            # generate a directory for each fragment
            Path(f"fragment_{i}").mkdir(exist_ok=True)
            fragmol.write_xyz_to_file(f"fragment_{i}/fragment_{i}.xyz")

    # if first fragment has less atoms than min_nat or more than max_nat, raise an error
    # NOTE: MM Check only for min_nat,
    # as fragment can still fragment into smaller fragments at the next step
    if fragmols[0].num_atoms < min_nat:
        raise RuntimeError(
            "First fragment has less atoms than the minimum number of atoms."
        )
    # Optimize the first fragment
    try:
        optfragmol = engine.optimize(fragmols[0])
    except RuntimeError as e:
        raise RuntimeError(f"Fragment optimization failed: {e}") from e
    # Differentiate again in fragments and take only the first one
    optfragmol = detect_fragments(optfragmol, verbosity)[0]

    # Check again if the fragment has less atoms than min_nat or more than max_nat, raise an error
    if optfragmol.num_atoms < min_nat or optfragmol.num_atoms > max_nat:
        raise RuntimeError("Fragment has less atoms than the minimum number of atoms.")

    # Check if the HL gap is larger than a given threshold
    if not engine.check_gap(optfragmol):
        raise RuntimeError("HL gap is smaller than the threshold.")

    return optfragmol


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
            sum_radii = get_cov_radii(mol.ati[i], COV_RADII) + get_cov_radii(
                mol.ati[j], COV_RADII
            )
            if distance <= sum_radii * 1:
                graph.add_edge(i, j)

    # Detect connected components (fragments)
    fragments = [list(component) for component in nx.connected_components(graph)]
    # sort the fragments by len
    fragments = sorted(fragments, key=len, reverse=True)

    # Generate a Molecule object for each fragment
    fragment_molecules = []
    for counter, fragment in enumerate(fragments):
        if verbosity > 1:
            print(f"Fragment: {counter}:\n{fragment}")
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
        fragment_molecule.charge = set_random_charge(fragment_molecule.ati, verbosity)
        if verbosity > 1:
            print(f"Fragment molecule: {fragment_molecule}")
        # Append the fragment molecule to the list
        fragment_molecules.append(fragment_molecule)

    return fragment_molecules


def get_cov_radii(at: int, vdw_radii: str = "mlmgen") -> float:
    """
    D3 covalent radii.
    """
    if vdw_radii == "mlmgen":
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
    elif vdw_radii == "pyykko":
        # Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
        # Values for metals decreased by 10%
        covalent_rad_2009 = [
            0.32,
            0.46,  # H, He
            1.20,
            0.94,
            0.77,
            0.75,
            0.71,
            0.63,
            0.64,
            0.67,  # Li-Ne
            1.40,
            1.25,
            1.13,
            1.04,
            1.10,
            1.02,
            0.99,
            0.96,  # Na-Ar
            1.76,
            1.54,  # K, Ca
            1.33,
            1.22,
            1.21,
            1.10,
            1.07,  # Sc-
            1.04,
            1.00,
            0.99,
            1.01,
            1.09,  # -Zn
            1.12,
            1.09,
            1.15,
            1.10,
            1.14,
            1.17,  # Ga-Kr
            1.89,
            1.67,  # Rb, Sr
            1.47,
            1.39,
            1.32,
            1.24,
            1.15,  # Y-
            1.13,
            1.13,
            1.08,
            1.15,
            1.23,  # -Cd
            1.28,
            1.26,
            1.26,
            1.23,
            1.32,
            1.31,  # In-Xe
            2.09,
            1.76,  # Cs, Ba
            1.62,
            1.47,
            1.58,
            1.57,
            1.56,
            1.55,
            1.51,  # La-Eu
            1.52,
            1.51,
            1.50,
            1.49,
            1.49,
            1.48,
            1.53,  # Gd-Yb
            1.46,
            1.37,
            1.31,
            1.23,
            1.18,  # Lu-
            1.16,
            1.11,
            1.12,
            1.13,
            1.32,  # -Hg
            1.30,
            1.30,
            1.36,
            1.31,
            1.38,
            1.42,  # Tl-Rn
            2.01,
            1.81,  # Fr, Ra
            1.67,
            1.58,
            1.52,
            1.53,
            1.54,
            1.55,
            1.49,  # Ac-Am
            1.49,
            1.51,
            1.51,
            1.48,
            1.50,
            1.56,
            1.58,  # Cm-No
            1.45,
            1.41,
            1.34,
            1.29,
            1.27,  # Lr-
            1.21,
            1.16,
            1.15,
            1.09,
            1.22,  # -Cn
            1.36,
            1.43,
            1.46,
            1.58,
            1.48,
            1.57,  # Nh-Og
        ]
        # D3 covalent radii used to construct the coordination number
        covalent_rad_d3 = [4.0 / 3.0 * rad for rad in covalent_rad_2009]
        return covalent_rad_d3[at]
    else:
        raise ValueError("Invalid vdw_radii argument.")
