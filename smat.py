from math import sqrt
import numpy
from Bio import pairwise2
from Bio.PDB import *
from Bio.PDB import Residue
from Bio.PDB import NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.SubsMat.MatrixInfo import blosum62
from pandas import read_csv, DataFrame
from pprint import pprint


# when computing electrostatic potential:
# - if target point is inside Van der Vaalse radius of current atom,
#   skip this atom

aa_residue_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
cofactor_residue_names = ['GDP', 'GTP', 'ADP', 'ATP', 'FMN', 'FAD', 'NAD', 'HEM']
aa_residue_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                      'V']
aa_residues = dict()


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


class ResidueDesc:
    def __init__(self, long_name: str, short_name: str, atoms=list()):
        self.long_name = long_name
        self.short_name = short_name
        self.atoms = atoms

    def __repr__(self):
        return f'Residue name: {self.long_name}\n' + \
            f'Residue short name: {self.short_name}\n' + \
            f'Residue atoms: {[x.get_type() for x in self.atoms]}'

    def get_long_name(self):
        return self.long_name

    def get_short_name(self):
        return self.short_name

    def get_atoms(self):
        return self.atoms


class AtomDesc:
    def __init__(self, type: str = '', charge = float(), edep = float(), A = float(), B = float(), mass = float(), parent = ResidueDesc):
        self.type = type
        self.charge = charge
        self.edep = edep
        self.A = A
        self.B = B
        self.mass = mass
        self.parent = parent

    def __repr__(self):
        return f'Atom type: {self.type}\n' + \
            f'Atom charge: {self.charge}\n' + \
            f'Atom edep: {self.edep}\n' + \
            f'Atom Van der Waals 12-part (A): {self.A}\n' + \
            f'Atom Van der Waals 6-part (B): {self.B}\n' + \
            f'Atom mass: {self.mass}\n' + \
            f'Atom parent: {self.parent.get_long_name()}'

    def get_type(self):
        return self.type

    def get_charge(self):
        return self.charge

    def get_edep(self):
        return self.edep

    def get_A(self):
        return self.A

    def get_B(self):
        return self.B

    def get_mass(self):
        return self.mass

    def get_parent(self):
        return self.parent


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


def point_belongs_to_active_site(point: list, atoms: list, center: list) -> bool:
    # for atom in atoms:
    #     radius = get_van_der_waals_radius(atom)
    #     if ray_intersects_sphere(center, point, atom.coord, radius):
    #         return False

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


def get_atoms_description() ->dict:
    # get atom types description ones for all other functions if needed
    # some aa have 2 variants of protonation

    # divide prep file on described residues per lines
    elements = dict()
    index = 0
    elements[index] = []
    with open('Docking_killer/database/prep/all.in', 'r') as file:
        for line in file.readlines()[2:]:
            if line.rstrip() != 'DONE' and line.rstrip() != 'STOP':
               elements[index].append(line.rstrip())
            elif line.rstrip() == 'DONE':
                index += 1
                elements[index] = []
                continue
            else:
                break

    # get some data from csv table
    data = read_csv('Docking_killer/VanDerWaals_2.csv', header=0, delimiter=';')

    # construct the final dict with proper data
    # parsing the dictionary with lines
    residues = list()
    for key in list(elements.keys())[:-1]:

        # find molecules with separated charges notes
        separate_charges = list() # if we have such - they will store here
        for line in elements[list(elements.keys())[key]]:
            if line.rstrip().split(' ')[0] == 'CHARGE':

                # get these fucking charges as array of strings
                index = elements[list(elements.keys())[key]].index(line) + 1
                while line != '':
                    line = elements[list(elements.keys())[key]][index]
                    separate_charges = separate_charges + [x for x in line.rstrip().split(' ') if x != '']
                    index += 1

        # create ResidueDesc object and fill the variables
        res_desc = ResidueDesc(long_name=elements[key][0], short_name=elements[key][2].split(' ')[1])
        for line in range(5, len(elements[key])):
            if elements[key][line] != '':
                line_elements = [x for x in elements[key][line].split(' ') if x != '']
                # filter the dummy atoms
                if line_elements and line_elements[1] != 'DUMM':

                    # filter the lone pairs
                    if line_elements[2] == 'LP':
                        continue

                    # create AtomDesc object for all atoms in residue
                    atom_desc = AtomDesc(type=line_elements[2])

                    # fill the charges
                    if len(separate_charges) == 0:
                        atom_desc.charge = float(line_elements[10])
                    else:
                        # check if dummy atoms have charges
                        #calculate the number of atoms
                        index = 0
                        for i in range(5, len(elements[key])):
                            if elements[key][i] != '':
                                index += 1
                            else:
                                break

                        if len(separate_charges) == index: # they have
                            atom_desc.charge = float(separate_charges[line - 5])
                        else: # they don't
                            atom_desc.charge = float(separate_charges[line - 8])

                    atom_type = ''
                    for atom in range(0, len(data)):
                        # fill the properties for common atoms
                        if data['Atom'][atom] == atom_desc.get_type():
                            atom_type = list(data['Atom']).index(atom_desc.get_type())
                        elif atom_desc.get_type()[0] == 'C' and atom_desc.get_type() != 'CZ':
                            atom_type = list(data['Atom']).index('C*')
                        elif atom_desc.get_type()[0] == 'N':
                            atom_type = list(data['Atom']).index('N')
                        elif atom_desc.get_type() == 'HW':
                            atom_type = list(data['Atom']).index('HO')

                    atom_desc.radius = float(data['Radius'][atom_type])
                    atom_desc.edep = float(data['Edep'][atom_type])
                    atom_desc.A = float(data['A'][atom_type])
                    atom_desc.B = float(data['B'][atom_type])
                    atom_desc.mass = float(data['Mass'][atom_type])

                    # fill the properties for atoms that are not normal
                    # we have undefined atom type CZ - for sp hybridisated C

                    res_desc.atoms.append(atom_desc)

            else:
                break

        # todo all residues have atoms from HEM - Why?!

        residues.append(res_desc)

    # for residue in residues:
    #     pprint(residue)

    return residues


def calculate_potential(atoms: list, atoms_desc: DataFrame) -> float:
    # point is a list of coordinates like [x, y, z]
    # atoms is a list of Atom objects

    print(atoms[123].get_parent())


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
    atoms_desc = get_atoms_description()
    # calculate_potential(neighbour_atoms, atoms_desc)


    class NeighbourSelect(Select):
        def accept_atom(self, atom):
            if atom in neighbour_atoms:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(structure)
    io.save('neighbours_only.pdb', NeighbourSelect())

