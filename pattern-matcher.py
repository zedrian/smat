import warnings
from functools import reduce
import math
import numpy
import operator
import time

from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pandas import read_csv
from Bio.PDB import Atom, PDBParser, Superimposer


# class Atom:
#     def __init__(self, name: str, x: float, y: float, z: float):
#         self.name = name
#         self.coords = [x, y, z]
#         self.element = name[0]
#
#     def __eq__(self, other):
#         return self.name == other.name and self.coords == other.coords
#
#     def __lt__(self, other):
#         return self.element < other.element or self.name < other.name or \
#                self.coords[0] < other.coords[0] or self.coords[1] < other.coords[1] or \
#                self.coords[2] < other.coords[2]
#
#     def __hash__(self):
#         return self.element.__hash__() + self.name.__hash__() + \
#                self.coords[0].__hash__() + self.coords[1].__hash__() + self.coords[2].__hash__()
#
#     def distance_to(self, other_atom) -> float:
#         squared_distance = (self.coords[0] - other_atom.coords[0]) ** 2 + \
#                            (self.coords[1] - other_atom.coords[1]) ** 2 + \
#                            (self.coords[2] - other_atom.coords[2]) ** 2
#         return math.sqrt(squared_distance)
#
#     def __repr__(self):
#         return f'{self.name}:[{self.coords[0]: 0.3f} {self.coords[1]: 0.3f} {self.coords[2]: 0.3f}]'


class Pattern:
    def __init__(self, pdb: str, name: str, sequence: str, atoms: list):
        self.pdb = pdb
        self.name = name
        self.sequence = sequence
        # ugly code to remove duplicates from atoms
        self.atoms = []
        for new_atom in atoms:
            already_added = False
            for atom in self.atoms:
                if all(atom.coord == new_atom.coord):
                    already_added = True
                    break
            if not already_added:
                self.atoms.append(new_atom)

    def __repr__(self):
        return 'Pattern { ' + \
               f'pdb: {self.pdb}, ' + \
               f'name: {self.name}, ' + \
               f'sequence: {self.sequence}, ' + \
               f'atoms: {self.atoms} ' + \
               '}'

    def sort_atoms(self, element_process_priority: list):
        self.atoms = sorted(self.atoms, key=lambda a: element_process_priority.index(a.get_id()[0]))


class Match:
    def __init__(self, matched_pattern_atoms: list, matched_inspecting_atoms: list):
        self.matched_pattern_atoms = matched_pattern_atoms
        self.matched_inspecting_atoms = matched_inspecting_atoms
        self.calculate_average_distance_error()

    def generate_plus_atom_matches(self, next_pattern_atom: Atom, inspecting_atoms: list,
                                   distance_error_coefficient: float) -> list:
        pattern_atoms = self.matched_pattern_atoms + [next_pattern_atom]
        pattern_atoms_distances = [numpy.linalg.norm(atom.coord - next_pattern_atom.coord) for atom in self.matched_pattern_atoms]

        matches = list()
        for atom in inspecting_atoms:
            # skip atom that is already presented in matched atoms
            if atom in self.matched_inspecting_atoms:
                continue

            # skip non-compatible atoms
            if atom.get_id()[0] != next_pattern_atom.get_id()[0]:
                continue

            # skip atoms if they do not fit by distance
            all_distances_are_acceptable = True
            for atom_index in range(len(pattern_atoms_distances)):
                acceptable_distance = pattern_atoms_distances[atom_index]
                if not distance_is_acceptable(self.matched_inspecting_atoms[atom_index], atom, acceptable_distance, distance_error_coefficient):
                    all_distances_are_acceptable = False
                    break

            if all_distances_are_acceptable:
                matches.append(Match(pattern_atoms, self.matched_inspecting_atoms + [atom]))
        return matches

    def calculate_average_distance_error(self):
        # superimposer = Superimposer()
        # superimposer.set_atoms(self.matched_pattern_atoms, self.matched_inspecting_atoms)
        # superimposer.rms
        distance_errors = list()
        for first_atom_index in range(len(self.matched_pattern_atoms)):
            for second_atom_index in range(first_atom_index+1, len(self.matched_pattern_atoms)):
                pattern_distance = numpy.linalg.norm(self.matched_pattern_atoms[first_atom_index] - self.matched_pattern_atoms[second_atom_index])
                inspecting_distance = numpy.linalg.norm(self.matched_inspecting_atoms[first_atom_index] - self.matched_inspecting_atoms[second_atom_index])
                distance_errors.append(math.fabs(inspecting_distance / pattern_distance - 1.0))
        self.average_distance_error = reduce(operator.add, distance_errors, 0) / len(self.matched_pattern_atoms)

    def __repr__(self):
        return 'Match { ' + \
            f'matched_pattern_atoms: {[f"{a.get_id()}:{a.coord}" for a in self.matched_pattern_atoms]}, ' + \
            f'matched_inspecting_atoms: {[f"{a.get_id()}[{a.get_parent().get_id()[1]}]:{a.coord}" for a in self.matched_inspecting_atoms]}, ' + \
            f'average_distance_error: {self.average_distance_error} ' + \
            '}'


def distance_is_acceptable(atom_a: Atom, atom_b: Atom, acceptable_distance: float,
                           distance_error_coefficient: float) -> bool:
    distance = numpy.linalg.norm(atom_a.coord - atom_b.coord)
    distance_fraction = distance / acceptable_distance

    return 1 - distance_error_coefficient <= distance_fraction <= 1 + distance_error_coefficient


# returns list of Bio.PDB.Atom objects
def get_inspecting_atoms(path: str, element_process_priority: list) -> list:
    structure = PDBParser(QUIET=True).get_structure('structure', 'Data/CHO_surface.pdb')

    # get Bio.PDB.Atom objects
    atoms = [a for a in structure.get_atoms()]
    # remove unknown atoms
    atoms = [a for a in atoms if a.get_id()[0] in element_process_priority]
    # sort atoms by element process priority
    atoms = sorted(atoms, key=lambda a: element_process_priority.index(a.get_id()[0]))

    return atoms


# returns list of Atom objects
def get_pattern_atoms(atoms_string: str, element_process_priority: list) -> list:
    if atoms_string == 'FILE_NOT_FOUND':
        return list()

    atoms = list()

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)
        for atom_string in atoms_string.split(', '):
            components = atom_string.split(':')
            name = components[0]
            coords_string = components[1][1:-1]
            coords = [float(x) for x in coords_string.split(' ') if x != '']
            atoms.append(Atom.Atom(name, numpy.asarray(coords), None, None, name, name, None))
    atoms = sorted(atoms, key=lambda a: element_process_priority.index(a.get_id()[0]))

    return atoms


def create_pattern_from_row(row, element_process_priority: list) -> Pattern:
    pdb = row['PDB_code']
    name = row['Pattern_name']
    sequence = row['Pattern_sequence']
    atoms = get_pattern_atoms(row['Contacts'], element_process_priority)

    return Pattern(pdb=pdb, name=name, sequence=sequence, atoms=atoms)


# returns list of Match objects
def find_pattern_matches(pattern: Pattern, inspecting_atoms: list, element_process_priority: list,
                         distance_error_coefficient: float = 0.1) -> list:
    matches = list()

    pattern_atom_count = len(pattern.atoms)
    print(f'atoms in pattern: {pattern_atom_count}')

    if pattern_atom_count >= 1:
        first_pattern_atom = pattern.atoms[0]
        matches = [Match([first_pattern_atom], [a]) for a in inspecting_atoms if
                            a.get_id()[0] == pattern.atoms[0].get_id()[0]]
        print(f'1-atom matches: {len(matches)}')

        for atom_index in range(1, len(pattern.atoms)):
            next_pattern_atom = pattern.atoms[atom_index]
            next_matches = reduce(operator.add, [
                match.generate_plus_atom_matches(next_pattern_atom, inspecting_atoms, distance_error_coefficient)
                for match in matches
            ], [])
            print(f'{atom_index+1}-atoms matches: {len(next_matches)}')
            matches = next_matches

    return matches


if __name__ == '__main__':
    start_time = time.time()

    element_process_priority = ['S', 'N', 'O', 'C', 'H', 'P', 'F', 'B']
    patterns = read_csv('pdb-hits.extended.5.csv')

    inspecting_atoms = get_inspecting_atoms('Data/CHO_surface.pdb', element_process_priority)
    total_matches = []

    for row_index, pattern_row in patterns.iterrows():
        pattern = create_pattern_from_row(pattern_row, element_process_priority)
        if len(pattern.atoms) == 0:
            total_matches.append([])
            continue
        print(f'processing pattern {pattern.name} from {pattern.pdb}')
        matches = find_pattern_matches(pattern, inspecting_atoms, element_process_priority)
        total_matches.append(matches)

        print(f'{row_index+1} / {patterns.shape[0]}')

    patterns['Matches'] = total_matches
    patterns.to_csv('matches.csv')

    print(f'done in {time.time() - start_time:0.3} seconds')
