from math import sqrt
import numpy
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

van_der_waals_radiuses = {'H': 0.6000, 'HO': 0.0000, 'HS': 0.6000, 'HC': 1.4870, 'H1': 1.3870, 'H2': 1.2870, 'H3': 1.1870,
                          'HP': 1.1000, 'HA': 1.4590, 'H4': 1.4090, 'H5': 1.3590, 'HW': 0.0000, 'HZ': 1.4590, 'O': 1.6612,
                          'O2': 1.6612, 'OW': 1.7683, 'OH': 1.7210, 'OS': 1.6837, 'OP': 1.8500, 'C*': 1.9080, 'CI': 1.9080,
                          'C5': 1.9080, 'C4': 1.9080, 'CT': 1.9080, 'CX': 1.9080, 'C': 1.9080, 'N': 1.8240, 'N3': 1.8240,
                          'S': 2.0000, 'SH': 2.0000, 'P': 2.1000, 'MG': 0.7926, 'C0': 1.7131, 'Zn': 1.10, 'F': 1.75, 'Cl': 1.948,
                          'Br': 2.22, 'I': 2.35, 'EP': 0.00}


class BoundingBox:
    def __init__(self):
        self.min_x = 1.0e+10
        self.max_x = -1.0e+10
        self.min_y = 1.0e+10
        self.max_y = -1.0e+10
        self.min_z = 1.0e+10
        self.max_z = -1.0e+10

    def store_point(self, x: float, y: float, z: float):
        self.min_x = min(self.min_x, x)
        self.max_x = max(self.max_x, x)
        self.min_y = min(self.min_y, y)
        self.max_y = max(self.max_y, y)
        self.min_z = min(self.min_z, z)
        self.max_z = max(self.max_z, z)

    def __str__(self):
        return f'[{self.min_x}..{self.max_x}] [{self.min_y}..{self.max_y}] [{self.min_z}..{self.max_z}]'


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


def get_center_of_mass(residue: Residue) -> list:
    center_of_mass = None
    mass = 0.0
    for atom in residue.get_atoms():
        if center_of_mass is None:
            center_of_mass = atom.coord * atom.mass
        else:
            center_of_mass = center_of_mass + atom.coord * atom.mass
        mass = mass + atom.mass
    center_of_mass = center_of_mass / mass

    return [center_of_mass[0], center_of_mass[1], center_of_mass[2]]


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
    chain_center_of_mass = [chain_center_of_mass[0], chain_center_of_mass[1], chain_center_of_mass[2]]

    closest_ligand = None
    closest_ligand_to_center_of_mass_squared_distance = 1e30 # just a big number

    # find closest ligand
    for residue in chain.get_residues():
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
    print(f'ligand selected: {closest_ligand.get_resname()} (distance to chain CoM: {sqrt(closest_ligand_to_center_of_mass_squared_distance)})')

    return closest_ligand


def get_neighbor_atoms(chain: Chain, ligand: Residue) -> list:
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


def get_bounding_box(atoms: list) -> BoundingBox:
    box = BoundingBox()

    for atom in atoms:
        box.store_point(atom.coord[0], atom.coord[1], atom.coord[2])

    return box


def get_van_der_walls_radius(atom: Atom) -> float:
    if atom.get_id() in van_der_waals_radiuses:
        print(f'attempt to get van der Waals radius for unexpected atom: {atom.get_id()}')
        exit(1)
    return van_der_waals_radiuses[atom.get_id()]


def point_belongs_to_active_site(point: list, atoms: list, center: list) -> bool:
    def ray_intersects_sphere(ray_start: list, ray_end: list, sphere_center: list, sphere_radius: float) -> bool:
        return False

    max_distance_from_center = 20.0
    
    # is point too far from center?
    distance = sqrt((point[0] - center[0])**2 +
                    (point[1] - center[1])**2 +
                    (point[2] - center[2])**2)
    if distance > max_distance_from_center:
        return False

    # check whether center-to-point ray intersects any van der Waals radius of atoms
    for atom in atoms:
        radius = get_van_der_walls_radius(atom)
        if ray_intersects_sphere(center, point, atom.coord, radius):
            return False

    return True


def get_potential_grid_coordinates(neighbour_atoms: list, bounding_box: BoundingBox, ligand_center_of_mass: list) -> list:
    step = 0.2

    grid_coordinates = list()

    x = bounding_box.min_x
    while x < bounding_box.max_x:
        y = bounding_box.min_y
        while y < bounding_box.max_y:
            z = bounding_box.min_z
            while z < bounding_box.max_z:
                # do stuff
                current_point = [x, y, z]
                if point_belongs_to_active_site(current_point, neighbour_atoms, ligand_center_of_mass):
                    grid_coordinates.append(current_point)

                z += step

            y += step

        x += step

    return grid_coordinates


def calculate_potential(point: list, atoms: list) -> float:
    # point is a list of coordinates like [x, y, z]
    # atoms is a list of Atom objects

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
    neighbour_atoms = get_neighbor_atoms(chain, ligand)
    bounding_box = get_bounding_box(neighbour_atoms)
    grid_coordinates = get_potential_grid_coordinates(neighbour_atoms, bounding_box, get_center_of_mass(ligand))

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