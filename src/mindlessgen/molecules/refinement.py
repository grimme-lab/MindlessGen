"""
This module handles all optimization and fragment detection steps
to obtain finally a valid molecule.
"""

from threading import Event
import warnings
from pathlib import Path
import networkx as nx  # type: ignore
import numpy as np

from ..qm.base import QMMethod
from ..prog import GenerateConfig, RefineConfig, ResourceMonitor
from .molecule import Molecule
from .miscellaneous import (
    get_cov_radii,
    set_random_charge,
    calculate_protons,
    calculate_ligand_electrons,
    calculate_uhf,
    get_lanthanides,
    get_actinides,
)

COV_RADII = "pyykko"


def iterative_optimization(
    mol: Molecule,
    engine: QMMethod,
    config_generate: GenerateConfig,
    config_refine: RefineConfig,
    resources_local: ResourceMonitor,
    stop_event: Event,
    verbosity: int = 1,
) -> Molecule | None:
    """
    Iterative optimization and fragment detection.
    """
    rev_mol = mol.copy()
    previous_fragments = None  # To store atom counts from the previous cycle

    if config_refine.debug:
        verbosity = 3

    for cycle in range(config_refine.max_frag_cycles):
        # Run single points first, start optimization if scf converges
        try:
            with resources_local.occupy_cores(1):
                if stop_event.is_set():
                    return None
                _ = engine.singlepoint(rev_mol, 1, verbosity)
        except RuntimeError as e:
            raise RuntimeError(
                f"Single-point calculation failed at fragmentation cycle {cycle}: {e}"
            ) from e

        # Optimize the current molecule
        try:
            with resources_local.occupy_cores(config_refine.ncores):
                if stop_event.is_set():
                    return None
                rev_mol = engine.optimize(
                    rev_mol, config_refine.ncores, None, verbosity
                )
        except RuntimeError as e:
            raise RuntimeError(
                f"Optimization failed at fragmentation cycle {cycle}: {e}"
            ) from e

        if verbosity > 2:
            # Print coordinates of optimized molecule
            print(f"Optimized molecule in cycle {cycle + 1}:\n{rev_mol.xyz}")

        # Detect fragments from the optimized molecule
        fragmols = detect_fragments(
            mol=rev_mol,
            molecular_charge=config_generate.molecular_charge,
            vdw_scaling=config_generate.scale_fragment_detection,
            verbosity=verbosity,
        )

        # Extract the number of atoms from each fragment for comparison
        current_atom_counts = [fragmol.num_atoms for fragmol in fragmols]

        # Check if the current atom counts are the same as in the previous cycle
        # or if there is only one fragment
        if previous_fragments is not None:
            if current_atom_counts == previous_fragments:
                if verbosity > 0:
                    print(
                        "Fragments are identical to the previous cycle by atom count. "
                        + f"Stopping at cycle {cycle + 1}."
                    )
                break  # Stop if the atom counts are the same as in the previous cycle
        if len(fragmols) == 1:
            if verbosity > 1:
                print(
                    f"Only one fragment detected in cycle {cycle + 1}. "
                    + f"Stopping at cycle {cycle + 1}."
                )
            break
        if verbosity > 1:
            print(f"Fragmentation cycle {cycle + 1}: Fragments: {current_atom_counts}")

        previous_fragments = current_atom_counts  # Update for the next cycle

        if verbosity > 2:
            for i, fragmol in enumerate(fragmols):
                # Generate a directory for each fragment
                Path(f"cycle_{cycle + 1}_fragment_{i}").mkdir(exist_ok=True)
                fragmol.write_xyz_to_file(
                    f"cycle_{cycle + 1}_fragment_{i}/fragment_{i}.xyz"
                )

        # Select the first fragment for the next cycle
        if fragmols[0].num_atoms < config_generate.min_num_atoms:
            raise RuntimeError(
                f"Largest fragment in cycle {cycle + 1} has fewer atoms than the minimum required."
            )
        # Check if the composition of the largest fragment
        # is still in line with the cfg.element_composition
        for ati, elem_range in config_generate.element_composition.items():
            min_limit, max_limit = elem_range
            count = fragmols[0].atlist[ati]
            if min_limit is not None and count < min_limit:
                raise RuntimeError(
                    f"Element {ati} is underrepresented "
                    + f"in the largest fragment in cycle {cycle + 1}."
                )
            # Check is actually nonsense in the current workflow,
            # but we want to have this as agnostic as possible
            if max_limit is not None and count > max_limit:
                raise RuntimeError(
                    f"Element {ati} is overrepresented "
                    + f"in the largest fragment in cycle {cycle + 1}."
                )
        if config_generate.molecular_charge is not None:
            protons = calculate_protons(fragmols[0].atlist)
            nel = protons - config_generate.molecular_charge
            f_elem = any(
                count > 0 and (i in get_lanthanides() or i in get_actinides())
                for i, count in enumerate(fragmols[0].atlist)
            )
            if f_elem:
                ligand_electrons = calculate_ligand_electrons(fragmols[0].atlist, nel)
                if ligand_electrons % 2 != 0:
                    raise RuntimeError(
                        f"Number of electrons in the largest fragment in cycle {cycle + 1} is odd."
                    )
            elif nel % 2 != 0:
                raise RuntimeError(
                    f"Number of electrons in the largest fragment in cycle {cycle + 1} is odd."
                )
        rev_mol = fragmols[
            0
        ]  # Set current_mol to the first fragment for the next cycle

    # Check the final fragment for staying in the min/max limits after self-consistent fragmentation
    if (
        rev_mol.num_atoms < config_generate.min_num_atoms
        or rev_mol.num_atoms > config_generate.max_num_atoms
    ):
        raise RuntimeError(
            "Final fragment has an invalid number of atoms (less than min or more than max)."
        )

    try:
        with resources_local.occupy_cores(1):
            if stop_event.is_set():
                return None
            gap_sufficient = engine.check_gap(
                molecule=rev_mol,
                threshold=config_refine.hlgap,
                ncores=1,
                verbosity=verbosity,
            )
    except NotImplementedError:
        warnings.warn("HOMO-LUMO gap check not implemented with this engine.")
    except RuntimeError as e:
        raise RuntimeError("HOMO-LUMO gap could not be checked.") from e
    if not gap_sufficient:
        raise RuntimeError("HOMO-LUMO gap does not meet the lower threshold.")

    return rev_mol


def detect_fragments(
    mol: Molecule,
    molecular_charge: int,
    vdw_scaling: float,
    verbosity: int = 1,
) -> list[Molecule]:
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
            if verbosity > 2:
                print(f"Distance between atom {i} and {j}: {distance:6.3f}")
                print(
                    f"Covalent radii of atom {i} and {j}, "
                    + "and the effective threshold: "
                    + f"{get_cov_radii(mol.ati[i], COV_RADII):6.3f}, "
                    + f"{get_cov_radii(mol.ati[j], COV_RADII):6.3f}, "
                    + f"{(sum_radii * vdw_scaling):6.3f}"
                )
            if distance <= sum_radii * vdw_scaling:
                graph.add_edge(i, j)

    # Detect connected components (fragments)
    fragments = [list(component) for component in nx.connected_components(graph)]
    # sort the fragments by len
    fragments = sorted(fragments, key=len, reverse=True)

    # Generate a Molecule object for each fragment
    fragment_molecules = []
    # if there is only one fragment, return the molecule as is
    if len(fragments) == 1:
        fragment_molecule = mol.copy()
        fragment_molecules.append(fragment_molecule)
        return fragment_molecules
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
        fragment_molecule.atlist = np.zeros(103, dtype=int)
        for atom in fragment_molecule.ati:
            fragment_molecule.atlist[atom] += 1
        # Update the charge of the fragment molecule
        if molecular_charge is not None:
            fragment_molecule.charge = molecular_charge
            fragment_molecule.uhf = calculate_uhf(fragment_molecule.atlist)
        else:
            fragment_molecule.charge, fragment_molecule.uhf = set_random_charge(
                fragment_molecule.ati, verbosity
            )
        fragment_molecule.set_name_from_formula()
        if verbosity > 1:
            print(f"Fragment molecule: {fragment_molecule}")
        # Append the fragment molecule to the list
        fragment_molecules.append(fragment_molecule)

    return fragment_molecules
