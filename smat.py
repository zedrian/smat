from math import sqrt
from Bio import pairwise2
from Bio.PDB import *
from Bio.PDB import Residue
from Bio.PDB import NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.SubsMat.MatrixInfo import blosum62


# when computing electrostatic potential:
# - if target point is inside Van der Vaalse radius of current atom,
#   skip this atom

aa_residue_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
cofactor_residue_names = ['GDP', 'GTP', 'ADP', 'ATP', 'FMN', 'FAD', 'NAD', 'HEM']
aa_residue_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                      'V']
aa_residues = dict()


def get_chain(structure: Structure) -> Chain:
    chains = [c for c in structure.get_chains()]

    # fast return if structure contains a single chain
    if len(chains) == 1:
        print('structure contains a single chain')
        return chains[0]
    print(f'structure contains {len(chains)} chains')

    def chain_contains_ligands(chain: Chain) -> bool:
        for residue in chain.get_residues():
            if not residue.get_resname() in aa_residue_names and not residue.get_resname() in cofactor_residue_names:
                atoms = [a for a in residue.get_atoms() if not a.get_name().startswith('H')]
                if len(atoms) >= 6:
                    return True
        return False
    chains_with_ligands = list(filter(chain_contains_ligands, chains))
    print(f'chains with ligands: {len(chains_with_ligands)}')

    def get_shortened_sequence(chain: Chain) -> str:
        line = ''
        for residue in chain.get_residues():
            if not residue.get_resname() in aa_residues.keys():
                continue
            line += aa_residues[residue.get_resname()]
        return line
    shortened_chain_sequences = list(map(get_shortened_sequence, chains_with_ligands))

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
    def get_chain_length(chain: Chain) -> int:
        sequence = get_shortened_sequence(chain)
        return len(sequence)
    sorted_chains = sorted(chains_with_ligands, key=get_chain_length, reverse=True)

    longest_chain = sorted_chains[0]
    print(f'chain selected: {longest_chain.get_id()}')

    return longest_chain


def is_ligand(residue: Residue) -> bool:
    if not residue.get_resname() in aa_residue_names and not residue.get_resname() in cofactor_residue_names:
        atoms = [a for a in residue.get_atoms() if not a.get_name().startswith('H')]
        if len(atoms) >= 6:
            return True


def get_ligand(chain: Chain) -> Residue:
    # - go through residues that are not AA and cofactors and are greater than of 6 atoms (not including H)
    # - all of them are candidates for being ligands
    # - find center of masses for the whole chain
    # - get candidate that is nearest to the center

    # compute chain center of mass at first
    chain_center_of_mass = None
    chain_mass = 0.0
    for residue in chain.get_residues():
        for atom in residue.get_atoms():
            if chain_center_of_mass is None:
                chain_center_of_mass = atom.coord * atom.mass
            else:
                chain_center_of_mass = chain_center_of_mass + atom.coord * atom.mass
            chain_mass = chain_mass + atom.mass
    chain_center_of_mass = chain_center_of_mass / chain_mass

    closest_ligand = None
    closest_ligand_to_center_of_mass_squared_distance = 1e30 # just a big number

    # find closest ligand
    for residue in chain.get_residues():
        if not is_ligand(residue):
            continue

        # compute ligand's center of mass
        ligand_center_of_mass = None
        ligand_mass = 0.0
        for atom in residue.get_atoms():
            if ligand_center_of_mass is None:
                ligand_center_of_mass = atom.coord * atom.mass
            else:
                ligand_center_of_mass = ligand_center_of_mass + atom.coord * atom.mass
            ligand_mass = ligand_mass + atom.mass
        ligand_center_of_mass = ligand_center_of_mass / ligand_mass

        # measure distance to chain center of mass
        delta_position = ligand_center_of_mass - chain_center_of_mass
        squared_distance = delta_position[0]**2 + delta_position[1]**2 + delta_position[2]**2

        # compare with current best result
        if squared_distance < closest_ligand_to_center_of_mass_squared_distance:
            closest_ligand = residue
            closest_ligand_to_center_of_mass_squared_distance = squared_distance

    # show result
    print(f'ligand selected: {closest_ligand.get_resname()} (distance to chain CoM: {sqrt(closest_ligand_to_center_of_mass_squared_distance)})')

    return closest_ligand


def find_contacts(chain: Chain, ligand: Residue) -> list:
    # use biopython neighboursearch to get list of AA atoms that close enough to ligand's atoms (using 10 angstroms)
    # - get list of all chain's atoms except ligand's atoms
    # - for each ligand's atom run neighbor search to find neighbors
    # - add that neighbors to result list of neighbors (exclude duplicates)

    # collect chain atoms
    chain_atoms = list()
    for residue in chain.get_residues():
        # do not count atoms from ligand itself
        if residue.get_resname() == ligand.get_resname():
            continue

        for atom in residue.get_atoms():
            chain_atoms.append(atom)

    neighbour_atoms = list()
    for atom in ligand.get_atoms():
        search = NeighborSearch(chain_atoms)
        current_neighbours = search.search(atom.coord, 10.0)
        for neighbour in current_neighbours:
            if neighbour not in neighbour_atoms:
                neighbour_atoms.append(neighbour)
    return neighbour_atoms


def get_potential_grid_coordinates(neighbour_atoms: list) -> list:
    #
    pass


def fill_potential_grid():
    # https://github.com/choderalab/ambermini/blob/master/share/amber/dat/leap/parm/parm10.dat
    # http://ambermd.org/formats.html#parm.dat
    pass

if __name__ == '__main__':
    # prepare AA residues dictionary
    for i in range(len(aa_residue_names)):
        aa_residues[aa_residue_names[i]] = aa_residue_letters[i]

    parser = PDBParser()
    structure = parser.get_structure('6b82', 'Docking_killer/proteins/CYPs/6b82_referense.pdb')
    chain = get_chain(structure)
    ligand = get_ligand(chain)
    neighbour_atoms = find_contacts(chain, ligand)

    # with open('neighbours.csv', 'w') as file:
    #     file.write('x,y,z\n')
    #     for atom in neighbour_atoms:
    #         file.write(f'{atom.coord[0]},{atom.coord[1]},{atom.coord[2]}\n')

    class NeighbourSelect(Select):
        def accept_atom(self, atom):
            if atom in neighbour_atoms:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(structure)
    io.save('neighbours_only.pdb', NeighbourSelect())

    print(len(neighbour_atoms))
