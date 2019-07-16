from math import sqrt
import numpy
from Bio import pairwise2
from Bio.PDB import *
from Bio.PDB import Residue
from Bio.PDB import NeighborSearch
from Bio.PDB.Atom import Atom
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

van_der_waals_radiuses = {'Br': 2.22, 'C': 1.908, 'C*': 1.908, 'CA': 1.908, 'CB': 1.908, 'CC': 1.908, 'CD': 1.908,
                          'CI': 1.908, 'CK': 1.908, 'CP': 1.908, 'CM': 1.908, 'CS': 1.908, 'CN': 1.908, 'CQ': 1.908,
                          'CR': 1.908, 'CV': 1.908, 'CW': 1.908, 'CY': 1.908, 'C0': 1.7131, 'CZ': 1.908, 'C5': 1.908,
                          'C4': 1.908, 'CT': 1.908, 'CX': 1.908, 'Cl': 1.948, 'EP': 0.0, 'F': 1.75, 'I': 2.35, 'H': 0.6,
                          'HO': 0.0, 'HS': 0.6, 'HC': 1.487, 'H1': 1.387, 'H2': 1.287, 'H3': 1.187, 'HP': 1.1,
                          'HA': 1.459, 'H4': 1.409, 'H5': 1.359, 'HW': 0.0, 'HZ': 1.459,'MG': 0.7926, 'N': 1.824,
                          'N3': 1.84, 'O': 1.6612, 'O2': 1.6612, 'OD': 1.6612, 'OW': 1.7683, 'OH': 1.721, 'OS': 1.6837,
                          'OP': 1.85, 'P': 2.1, 'S': 2.0, 'SH': 2.0,  'Zn': 1.1,
                          # next ones are untrusted:
                          'CE': 1.908, 'CG': 1.908, 'NE': 1.824, 'NH': 1.824, 'OE': 1.6612, 'OG': 1.6612, 'SD': 2.0,
                          'SG': 2.0, 'OXT': 1.7683}


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


class AtomDesc:
    def __init__(self, fullname: str = '', type: str = '', charge = float(), edep = float(), A = float(), B = float(), mass = float(), parent_name: str = ''):
        self.fullname = fullname
        self.type = type
        self.charge = charge
        self.edep = edep
        self.A = A
        self.B = B
        self.mass = mass
        self.parent_name = parent_name

    def __repr__(self):
        return f'Atom full name: {self.fullname}\n' + \
            f'Atom type: {self.type}\n' + \
            f'Atom charge: {self.charge}\n' + \
            f'Atom edep: {self.edep}\n' + \
            f'Atom Van der Waals 12-part (A): {self.A}\n' + \
            f'Atom Van der Waals 6-part (B): {self.B}\n' + \
            f'Atom mass: {self.mass}\n' + \
            f'Atom parent: {self.parent_name}'

    def get_fullname(self):
        return self.fullname

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

    def get_parent_name(self):
        return self.parent_name


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

    def get_atom(self, atom_fullname: str) -> AtomDesc:
        for atom in self.atoms:
            if atom.get_fullname() == atom_fullname:
                return atom

        print(f'ERROR: attempt to get atom {atom_fullname} from residue {self.long_name}')
        return None


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
    id = atom.get_id()

    if id not in van_der_waals_radiuses:
        if len(id) > 2:
            id = id[0:2]

    if id not in van_der_waals_radiuses:
        print(f'attempt to get van der Waals radius for unexpected atom: {atom.get_id()}')
        exit(1)
    return van_der_waals_radiuses[id]


def get_length(v: list) -> float:
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def point_belongs_to_active_site(point: list, atoms: list, center: list) -> bool:
    def get_direction(start: list, end: list) -> list:
        delta = numpy.subtract(end, start)
        length = get_length(delta)
        return [delta[0] / length, delta[1] / length, delta[2] / length]

    def ray_intersects_sphere(ray_start: list, ray_end: list, sphere_center: list, sphere_radius: float) -> bool:
        ray_length = get_length(numpy.subtract(ray_start, ray_end))

        # find sphere center projection on ray
        ray_direction = get_direction(ray_start, ray_end)
        time = numpy.dot(numpy.subtract(sphere_center, ray_start), ray_direction)
        sphere_center_projection = numpy.add(ray_start, [ray_direction[0]*time, ray_direction[1]*time,
                                                         ray_direction[2]*time])

        # measure distance from ray to sphere center:
        # - if greater than radius, return False
        sphere_center_to_ray_distance = get_length(numpy.subtract(sphere_center_projection, sphere_center))
        if sphere_center_to_ray_distance > sphere_radius:
            return False

        # check projection time:
        # - if between start and end, return True
        # - otherwise measure distance to start and end points
        if 0 < time < ray_length:
            return True
        if time < 0:
            return get_length(numpy.subtract(sphere_center, ray_start)) < sphere_radius
        return get_length(numpy.subtract(sphere_center, ray_end)) < sphere_radius

    max_distance_from_center = 20.0

    # is point too far from center?
    distance = get_length(numpy.subtract(point, center))
    if distance > max_distance_from_center:
        print('too far')
        return False

    # check whether center-to-point ray intersects any van der Waals radius of atoms
    for atom in atoms:
        radius = get_van_der_walls_radius(atom)
        if ray_intersects_sphere(center, point, atom.coord, radius):
            print('intersects')
            return False

    print('don\'t intersects')
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
    data = read_csv('Docking_killer/VanDerWaals.csv', header=0, delimiter=';')

    # construct the final dict with proper data
    # parsing the dictionary with lines
    residues = dict()
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
        res_desc = ResidueDesc(long_name=elements[key][0], short_name=elements[key][2].split(' ')[1], atoms=list())
        for line in range(5, len(elements[key])):
            if elements[key][line] != '':
                line_elements = [x for x in elements[key][line].split(' ') if x != '']
                # filter the dummy atoms
                if line_elements and line_elements[1] != 'DUMM':

                    # filter the lone pairs
                    if line_elements[2] == 'LP':
                        continue

                    # create AtomDesc object for all atoms in residue
                    atom_desc = AtomDesc(fullname=line_elements[1], type=line_elements[2], parent_name=res_desc.get_short_name())

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

        residues[res_desc.get_short_name()] = res_desc

    return residues


def calculate_potential(point: list, atoms: list, residues: dict) -> (float, float):
    # point is a list of coordinates like [x, y, z]
    # atoms is a list of Atom objects
    # residues is a dictionary with residue descriptions

    total_coulomb_potential = 0.0
    total_lennard_jones_energy = 0.0

    for atom in atoms:
        residue_name = atom.get_parent().get_resname()
        atom_description = residues[residue_name].get_atom(atom.get_fullname())

        # calculate Coulomb potential
        charge = atom_description.get_charge()
        distance = get_length(numpy.subtract(point, atom.get_coord()))
        coulomb_potential = charge / distance**2
        total_coulomb_potential += coulomb_potential

        # check whether we should calculate Lennard-Jones potential
        sigma = (atom_description.get_radius() + 2.35) / 2  # 2.35 is VdW-radius for Iodine (the greatest possible one)
        if distance > sigma * 2.5:  # 2.5 sigma is critical distance
            continue

        epsilon = atom_description.get_edep()

        lennard_jones_energy = 4*epsilon*((sigma / distance)**12 - (sigma / distance)**6)
        total_lennard_jones_energy += lennard_jones_energy

    return (total_coulomb_potential, total_lennard_jones_energy)


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

