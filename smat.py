from database_parser import get_residues_description, database
from math import sqrt, isnan, pi
import numpy
from sys import stdout, argv
from Bio.PDB import *
from Bio.PDB.Atom import Atom
from chainDesc  import ChainDesc
from classes import BoundingBox, ResiduesDatabase
from pdb_parser import *
from files_manager import *


# smart progress bar
def show_progress(label, width, percentage):
    progress = '['
    for i in range(0, width):
        if i / width < percentage:
            progress += '#'
        else:
            progress += ' '
    progress += '] {0:.3%}'.format(percentage)
    print('\r' + label + progress, end='')
    if percentage == 1.0:
        print('')
    stdout.flush()


def get_length(v: list) -> float:
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def get_direction(start: list, end: list) -> list:
    delta = numpy.subtract(end, start)
    length = get_length(delta)
    return [delta[0] / length, delta[1] / length, delta[2] / length]


def point_belongs_to_active_site(point: list, atoms: list, center: list) -> bool:

    def ray_intersects_sphere(ray_start: list, ray_end: list, sphere_center: list, sphere_radius: float) -> bool:
        ray_length = get_length(numpy.subtract(ray_start, ray_end))

        # find sphere center projection on ray
        ray_direction = get_direction(ray_start, ray_end)
        time = numpy.dot(numpy.subtract(sphere_center, ray_start), ray_direction)
        sphere_center_projection = numpy.add(ray_start, [ray_direction[0]*time, ray_direction[1]*time,
                                                         ray_direction[2]*time])

        # measure distance from ray to sphere center:
        # - if greater than radius, return False
        sctrd = numpy.subtract(sphere_center_projection, sphere_center)
        sphere_center_to_ray_distance = get_length(sctrd)
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
        return False

    # check whether center-to-point ray intersects any van der Waals radius of atoms
    for atom in atoms:  # physical atom
        radius = atom.get_atom_desc().get_radius()
        if ray_intersects_sphere(center, point, atom.get_coords(), radius):
            return False

    return True


def get_potential_grid_coordinates(step: float, neighbour_atoms: list, bounding_box: BoundingBox, ligand_center_of_mass: list, ligand_atoms: list) -> list:
    grid_coordinates = list()

    total_point_count = ((bounding_box.max_x - bounding_box.min_x)//step + 1) * \
                        ((bounding_box.max_y - bounding_box.min_y)//step + 1) * \
                        ((bounding_box.max_z - bounding_box.min_z)//step + 1)

    point_index = 0
    show_progress('potential grid calculation: ', 80, float(point_index)/float(total_point_count))
    x = bounding_box.min_x
    while x < bounding_box.max_x:
        column = list()
        y = bounding_box.min_y
        while y < bounding_box.max_y:
            row = list()
            z = bounding_box.min_z
            while z < bounding_box.max_z:
                # do stuff
                current_point = [x, y, z]
                if point_belongs_to_active_site(current_point, neighbour_atoms, ligand_center_of_mass):
                    row.append(current_point)
                else:
                    row.append(None)

                point_index += 1
                show_progress('potential grid calculation: ', 80, float(point_index) / float(total_point_count))

                z += step

            column.append(row)
            y += step

        grid_coordinates.append(column)
        x += step

    for atom in ligand_atoms:  # physical atoms
        coords = atom.get_coords()
        for column_index in range(len(grid_coordinates)):
            column = grid_coordinates[column_index]
            for row_index in range(len(column)):
                row = column[row_index]
                for point_index in range(len(row)):
                    point = row[point_index]
                    if point is not None:
                        if get_length(numpy.subtract(point, coords)) <= atom.get_atom_desc().get_radius():
                            grid_coordinates[column_index][row_index][point_index] = None

    show_progress('potential grid calculation: ', 80, 1.0)
    return grid_coordinates


def calculate_potential(point: list, atoms: list) -> (float, float):
    # point is a list of coordinates like [x, y, z]
    # atoms is a list of Atom objects

    total_coulomb_potential = 0.0
    total_lennard_jones_energy = 0.0

    atom_index = 0
    total_atom_count = len(atoms)
    show_progress('potential grid calculation: ', 80, float(atom_index) / float(total_atom_count))
    for atom in atoms:  # physical atom
        # calculate Coulomb potential
        charge = atom.get_atom_desc().get_charge()
        distance = get_length(numpy.subtract(point, atom.get_coords()))

        dielectric_const = 10
        coulomb_potential = charge / (4 * pi * dielectric_const * distance)
        total_coulomb_potential += coulomb_potential

        # check whether we should calculate Lennard-Jones potential
        sigma = atom.get_atom_desc().get_radius() * 2
        if distance <= sigma * 2.5:  # 2.5 sigma is critical distance
            epsilon = atom.get_atom_desc().get_edep()

            lennard_jones_energy = 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)

            total_lennard_jones_energy += lennard_jones_energy

        atom_index += 1
        show_progress('potential grid calculation: ', 80, float(atom_index) / float(total_atom_count))

    show_progress('potential grid calculation: ', 80, 1.0)

    return total_coulomb_potential, total_lennard_jones_energy


def construct_active_site_in_potentials_form(grid_coordinates: list, protein_atoms: list, ligand_atoms: list) -> list:
    # calculate indused potentials and Van Der Waals energy in each point of active center grid

    active_site_points = list()

    for column in grid_coordinates:
        for row in column:
            for point in row:
                if point is not None:
                    protein_coulomb_potential, protein_lennard_jones_energy = calculate_potential(point, protein_atoms)
                    ligand_coulomb_potential, ligand_lennard_jones_energy = calculate_potential(point, ligand_atoms)
                    if not isnan(protein_lennard_jones_energy):
                        if ligand_lennard_jones_energy != 0 and ligand_lennard_jones_energy < 10.0:
                            active_site_points.append(
                                {'coordinates': point, 'protein_coulomb': protein_coulomb_potential, f'protein_lennard_jones': protein_lennard_jones_energy,
                                'ligand_coulomb': ligand_coulomb_potential, 'ligand_lennard_jones': ligand_lennard_jones_energy})
                    else:
                        print(point)
    print('potentials calculated')

    return active_site_points


def calculate_forces(ligand_atoms: list, neighbour_atoms: list, dielectric_const: float = 10.0) -> dict:

    forces = dict()
    ligand_atom_count = float(len(ligand_atoms))
    progress = 0.0
    for ligand_atom in ligand_atoms:  # physical atom
        show_progress('calculating forces: ', 40, progress)
        integral_force = [0, 0, 0]
        integral_coulomb_force = [0, 0, 0]
        # TODO: refactor
        ligand_atom_charge = ligand_atom.get_atom_desc().get_charge()
        ligand_atom_edep = ligand_atom.get_atom_desc().get_edep()
        ligand_atom_radius = ligand_atom.get_atom_desc().get_radius()
        ligand_atom_position = ligand_atom.get_coords()
        for protein_atom in neighbour_atoms:  # physical atom
            protein_atom_charge = protein_atom.get_atom_desc().get_charge()
            protein_atom_edep = protein_atom.get_atom_desc().get_edep()
            protein_atom_radius = protein_atom.get_atom_desc().get_radius()
            protein_atom_position = protein_atom.get_coords()
            position_delta = numpy.subtract(protein_atom_position, ligand_atom_position)
            distance = get_length(position_delta)
            force_direction = get_direction(ligand_atom_position, protein_atom_position)

            if distance != 0.0 and distance < 10.1:
                # calculate coulomb force
                coulomb_force_module = (ligand_atom_charge * protein_atom_charge) / (4 * pi * dielectric_const * distance**2)

                lennard_force_module = 0.0
                if distance < ligand_atom_radius * 4.0:
                    # calculate the force according https://www.ks.uiuc.edu/Training/Workshop/SanFrancisco/lectures/Wednesday-ForceFields.pdf page 14
                    lennard_force_module = sqrt(ligand_atom_edep*protein_atom_edep)*(((ligand_atom_radius+protein_atom_radius)/distance)**12 - 2*((ligand_atom_radius+protein_atom_radius)/distance)**6)

                force_module = coulomb_force_module + lennard_force_module
                force = [force_direction[0] * force_module,
                                 force_direction[1] * force_module,
                                 force_direction[2] * force_module]
                integral_force = [integral_coulomb_force[0] + force[0],
                                  integral_coulomb_force[1] + force[1],
                                  integral_coulomb_force[2] + force[2]]

        # save to dictionary
        forces[ligand_atom.get_atom_desc().get_fullname()] = (ligand_atom_position, integral_force)

        progress += 1.0 / ligand_atom_count
    show_progress('calculating forces: ', 40, 1.0)

    return forces


def calculate_momentum(forces: dict, center_of_mass: list) -> list:
    total_momentum = [0.0, 0.0, 0.0]

    progress = 0.0
    show_progress('calculating momentum: ', 40, progress)
    for atom_name in forces:
        atom_position, force = forces[atom_name]
        radius = numpy.subtract(atom_position, center_of_mass)
        momentum = numpy.cross(radius, force)
        total_momentum = numpy.add(total_momentum, momentum)

        progress += 1.0 / len(forces)
        show_progress('calculating momentum: ', 40, progress)
    show_progress('calculating momentum: ', 40, 1.0)

    return total_momentum


def calculate_accumulated_force(forces: dict) -> list:
    accumulated_force = [0, 0, 0]
    for atom in forces:
        position, force = forces[atom]
        accumulated_force = [accumulated_force[0] + force[0],
                             accumulated_force[1] + force[1],
                             accumulated_force[2] + force[2]]
    return accumulated_force


if __name__ == '__main__':

    options, args = parse_input_command()
    if options.units == '':
        options.units = list()
    else:
        options.units = list(options.units)

    # get main variables from command line
    step = float(options.step)
    input_folder = options.folder

    # generate output folder
    results_folder = generate_output_folder(input_folder, step)

    # construct final dictionaries
    accumulated_forces = dict()
    momenta = dict()

    for root, dirs, filenames in os.walk(input_folder):
        file_count = float(len(filenames))
        file_index = 0.0
        for pdb_file in filenames:
            if not pdb_file.endswith('.pdb'):
                continue

            show_progress('processing PDB files: ', 40, file_index)
            print('')
            print(pdb_file)
            file_name = os.path.join(root, pdb_file)

            # get chain of the structure to work with
            parser = PDBParser()
            structure = parser.get_structure('pdb', file_name)
            chain = get_chain(structure)  # ChainDesc

            # get ligand and "transform" it to object of ResidueDesc class
            ligand = get_ligand(chain)  # Physical Residue
            ligand_atoms = ligand.get_atoms()  # Physical Atoms
            ligand_center_of_mass = get_center_of_mass(ligand)

            # write residues that need to be visualised to csv file
            units_short_names = ''.join(options.units).split(',')
            print(units_short_names)
            units = [u for u in chain.get_residues() if u.get_residue_desc().get_short_name() in units_short_names]  # physical residues
            units.append(ligand)  # add ligand (already a physical residue) to units
            write_units_csv(units, results_folder, pdb_file)

            neighbour_atoms = get_neighbor_atoms(chain, ligand)

            bounding_box = get_bounding_box(neighbour_atoms)

            forces = calculate_forces(ligand_atoms, neighbour_atoms)
            # save_forces_to_file(forces, results_folder, pdb_file)
            momentum = calculate_momentum(forces, ligand_center_of_mass)
            print(f'momentum = {momentum}')
            momenta[pdb_file] = momentum

            accumulated_force = calculate_accumulated_force(forces)
            # print(f'accumulated force = {forces}')
            accumulated_forces[pdb_file] = accumulated_force

            grid_coordinates = get_potential_grid_coordinates(step, neighbour_atoms, bounding_box, ligand_center_of_mass, ligand_atoms)
            print('grid coordinates calculated')

            active_site_points = construct_active_site_in_potentials_form(grid_coordinates, neighbour_atoms, ligand_atoms)
            save_potentials_to_file(active_site_points, results_folder, pdb_file)

            exit()
            file_index += 1.0 / file_count

    write_forces(results_folder, accumulated_forces)
    write_momenta(results_folder, momenta)

    class NeighbourSelect(Select):
        def accept_atom(self, atom):
            if atom in neighbour_atoms:
                return 1
            else:
                return 0