from Bio import pairwise2
from Bio.PDB import NeighborSearch
from Bio.PDB.Structure import Structure
from Bio.SubsMat.MatrixInfo import blosum62
from classes import BoundingBox, PhysicalResidue, PhysicalAtom
from chainDesc import ChainDesc
import numpy
from math import sqrt
from database_parser import database


def get_chain(structure: Structure) -> ChainDesc:
    chains = [ChainDesc(chain=c) for c in structure.get_chains()]

    # fast return if structure contains a single chain
    if len(chains) == 1:
        print('structure contains a single chain')
        return chains[0]
    print(f'structure contains {len(chains)} chains')

    chains_with_ligands = [c for c in chains if c.if_has_ligands()]
    print(f'chains with ligands: {len(chains_with_ligands)}')

    shortened_chain_sequences = [c.get_shortened_seq() for c in chains_with_ligands]

    # - check all pair alignments
    # - if at least one pair has equality score less than threshold,
    #   ask user for which chain to choose (by its letter)
    # - otherwise (that means that all chains are similar) choose longest one

    equality_threshold = 0.95

    for first_index in range(len(shortened_chain_sequences)):
        for second_index in range(first_index+1, len(shortened_chain_sequences)):
            first_sequence = shortened_chain_sequences[first_index]
            second_sequence = shortened_chain_sequences[second_index]
            alignments = pairwise2.align.globalds(first_sequence, second_sequence, blosum62, -10, -0.5)
            first_aligned, second_aligned, score, begin, end = alignments[0]
            if score < equality_threshold:
                print(f'two different sequences found (score={score})')
                print('please enter a letter of chain to work with: ', end='')
                chain_letter = input()[0]
                chain_index = int(chain_letter) - int('A')
                return chains[chain_index]

    # as we are here, then no different chains were found -
    # so choose longest one
    def get_chain_length(chain: ChainDesc) -> int:
        sequence = chain.get_shortened_seq()
        return len(sequence)

    sorted_chains = sorted(chains_with_ligands, key=get_chain_length, reverse=True)

    longest_chain = sorted_chains[0]
    print(f'chain selected: {longest_chain.chain.get_id()}')

    return longest_chain


def is_ligand(residue: PhysicalResidue) -> bool:
    if not residue.get_residue_desc().get_short_name() in database.get_amino_acids() + database.get_cofactors() or residue.get_residue_desc().check_deeper():
        atoms = [a.get_atom_desc() for a in residue.get_atoms() if not a.get_atom_desc().get_fullname().startswith('H')]
        if len(atoms) >= 6:
            return True


def get_center_of_mass(residue: PhysicalResidue) -> list:
    center_of_mass = None
    mass = 0.0
    for atom in residue.get_atoms():  # physical atom
        if center_of_mass is None:
            center_of_mass = atom.get_coords() * atom.get_atom_desc().get_mass()
        else:
            center_of_mass = center_of_mass + atom.get_coords() * atom.get_atom_desc().get_mass()
        mass = mass + atom.get_atom_desc().get_mass()
    center_of_mass = center_of_mass / mass

    return [center_of_mass[0], center_of_mass[1], center_of_mass[2]]


def get_ligand(chain: ChainDesc) -> PhysicalResidue:
    # - go through residues that are not AA and cofactors and are greater than of 6 atoms (not including H)
    # - all of them are candidates for being ligands
    # - find center of masses for the whole chain
    # - get candidate that is nearest to the center

    # compute chain center of mass at first
    chain_center_of_mass = None
    chain_mass = 0.0
    for residue in chain.get_residues():  # physical residue
        for atom in residue.get_atoms():  # physical atom
            if chain_center_of_mass is None:
                chain_center_of_mass = atom.coords * atom.get_atom_desc().get_mass()
            else:
                chain_center_of_mass = chain_center_of_mass + atom.coords * atom.get_atom_desc().get_mass()
            chain_mass = chain_mass + atom.get_atom_desc().get_mass()
    chain_center_of_mass = chain_center_of_mass / chain_mass
    chain_center_of_mass = [chain_center_of_mass[0], chain_center_of_mass[1], chain_center_of_mass[2]]

    closest_ligand = None
    closest_ligand_to_center_of_mass_squared_distance = 1e30  # just a big number

    # find closest ligand
    for residue in chain.get_residues():  # physical residue
        if not is_ligand(residue):
            continue

        ligand_center_of_mass = get_center_of_mass(residue)

        # measure distance to chain center of mass
        delta_position = numpy.subtract(ligand_center_of_mass, chain_center_of_mass)
        squared_distance = delta_position[0]**2 + delta_position[1]**2 + delta_position[2]**2

        # compare with current best result
        if squared_distance < closest_ligand_to_center_of_mass_squared_distance:
            closest_ligand = residue
            closest_ligand_to_center_of_mass_squared_distance = squared_distance

    # show result
    print(f'ligand selected: {closest_ligand.get_residue_desc().get_short_name()} distance to chain CoM: '
          f'{sqrt(closest_ligand_to_center_of_mass_squared_distance)})')

    return closest_ligand


def get_neighbor_atoms(chain: ChainDesc, ligand: PhysicalResidue) -> list:  # a list of physical atoms
    # use biopython neighboursearch to get list of AA atoms that close enough to ligand's atoms (using 10 angstroms)
    # - get list of all chain's atoms except ligand's atoms
    # - for each ligand's atom run neighbor search to find neighbors
    # - add that neighbors to result list of neighbors (exclude duplicates)

    # collect chain atoms
    chain_atoms = list()

    # get bio python residues of the chain
    bio_residues = list(chain.chain.get_residues())

    for residue in chain.get_residues():  # physical residue
        # do not count atoms from ligand itself
        if residue.get_residue_desc().get_short_name() == ligand.get_residue_desc().get_short_name():
            continue

        # do not count atoms from ligands
        # TODO: refactor
        if residue.get_residue_desc().get_short_name() not in database.get_amino_acids() + database.get_cofactors():
            continue
        elif residue.get_residue_desc().check_deeper():
            continue

        for atom in bio_residues[residue.get_index() - 1].get_atoms():  # BioPython atom!!!!
            chain_atoms.append(atom)

    neighbour_atoms = list()
    for atom in ligand.get_atoms():  # physical atom
        search = NeighborSearch(chain_atoms)
        current_neighbours = search.search(atom.get_coords(), 10.0)
        for neighbour in current_neighbours:  # BioPython atom!!!!
            if neighbour not in neighbour_atoms:
                neighbour_atoms.append(neighbour)

    neighbour_physical_atoms = list()
    for neighbour in neighbour_atoms:
        terminus = None
        if neighbour.get_parent().get_id()[1] == 1:
            terminus = 'N'
        elif 'OXT' in [a.get_id() for a in neighbour.get_parent().get_atoms()]:
            terminus = 'C'
        atom_desc = database.get_residue(neighbour.get_parent(), terminus).get_atom(neighbour.get_id())
        physical_atom = PhysicalAtom(bio_atom=neighbour, atom_desc=atom_desc, coords=neighbour.get_coord())
        neighbour_physical_atoms.append(physical_atom)

    return neighbour_physical_atoms


def get_bounding_box(atoms: list) -> BoundingBox:
    box = BoundingBox()

    for atom in atoms:  # physical atom
        box.store_point(atom.get_coords()[0], atom.get_coords()[1], atom.get_coords()[2])

    return box

