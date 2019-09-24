from optparse import OptionParser
from shutil import rmtree
import os
from smat import get_length


# parse the flags of the input command line
def parse_input_command():
    parser = OptionParser()
    parser.add_option("-i", "--file", dest="filename", help="input file to analyze (not realized yet)", metavar="FILE")  # todo make possible to run for 1 file
    parser.add_option("-f", "--folder", dest="folder", help="input folder with data", metavar="FOLDER")
    parser.add_option("-r", "--remaster", dest="remaster", default=False, help="remaster database (not realized yet)")  # todo generate database ones and write it to file and remaster each time it is necessary
    parser.add_option("-s", "--step", dest="step", help="the step of the grid", default=1.0, metavar="STEP")
    parser.add_option("-u", "--units", dest="units", help="residues to be saved", metavar="UNITS", default=None)
    parser.add_option("-k", "--dielectric_const", dest="diel_const", help="dielectric constant for inter-protein environment",  metavar="DIEL_CONST", default=10.0)

    return parser.parse_args()


# generate output folder and clear if already exists
def generate_output_folder(input_folder: str, step: float):
    results_folder = os.path.join(input_folder, f'results_{step}')
    if os.path.exists(results_folder):
        rmtree(results_folder)
    else:
        os.makedirs(results_folder)
        os.makedirs(os.path.join(results_folder, 'potentials'))
        os.makedirs(os.path.join(results_folder, 'units'))

    return results_folder


# write forces to file
def write_forces(output_folder: str, forces: dict):
    with open(os.path.join(output_folder, 'accumulated_forces.csv'), 'w') as file:
        file.write('pdb,force\n')
        for pdb in forces:
            file.write(f'{pdb},{get_length(forces[pdb])}\n')
        print('accumulated forces saved')


# write momenta to file
def write_momenta(output_folder: str, momenta: dict):
    with open(os.path.join(output_folder, 'momenta.csv'), 'w') as file:
        file.write('pdb,x,y,z,length\n')
        for pdb in momenta:
            momentum = momenta[pdb]
            file.write(f'{pdb},{momentum[0]},{momentum[1]},{momentum[2]},{get_length(momentum)}\n')
        print('momenta saved')


# write forces for each atom of ligand (atom name, coordinates, vectors)
def save_ligand_forces_to_file(forces: dict, results_folder: str, pdb_file: str):
    file_name = os.path.join(results_folder, f'{pdb_file[:-4]}_forces.csv')
    with open(file_name, 'w') as file:
        file.write('atom_name,x,y,z,force_x,force_y,force_z\n')
        for atom in forces:
            position, force = forces[atom]
            file.write(f'{atom},{position[0]},{position[1]},{position[2]},{force[0]},{force[1]},{force[2]}\n')
    print(f'forces saved to {file_name}')


# write units to csv file for the visualisator
def write_units_csv(units: list, results_folder: str, pdb_file: str):
    # create units folder if not exists
    units_folder = os.path.join(results_folder, 'units')
    if not os.path.exists(units_folder):
        os.makedirs(units_folder)

    with open(os.path.join(units_folder, f'{pdb_file[:-4]}.csv'), 'w') as file:
        file.write('x,y,z,radius,rgb,atom_name,atom_type,residue_name\n')
        for unit in units:  # physical residue
            if unit is not None:
                for atom in unit.get_atoms():  # physical atom
                    if atom.get_x() == atom.get_y() == atom.get_z() == 0.0:
                        continue
                    else:
                        rgb = '000000'
                        if atom.get_atom_desc().get_type()[0].upper() == 'C' and atom.get_atom_desc().get_type() != 'Cl':
                            rgb = 'FFFFCC'
                        elif atom.get_atom_desc().get_type()[0].upper() == 'N':
                            rgb = '3366FF'
                        elif atom.get_atom_desc().get_type()[0].upper() == 'O':
                            rgb = 'FF3300'
                        elif atom.get_atom_desc().get_type()[0].upper() == 'H':
                            rgb = 'FFCCFF'
                        elif atom.get_atom_desc().get_type()[0].upper() == 'S':
                            rgb = 'FFFF00'
                        elif atom.get_atom_desc().get_type() == 'Cl':
                            rgb = '33FF33'
                        elif atom.get_atom_desc().get_type() == 'Fe':
                            rgb = '996600'
                        line = f'{atom.get_x()},{atom.get_y()},{atom.get_z()},{atom.get_atom_desc().get_radius()/3},' +\
                            f'{rgb},{atom.get_atom_desc().get_fullname()},' + \
                            f'{atom.get_atom_desc().get_type()},{atom.get_atom_desc().get_parent_name()}\n'
                        file.write(line)
        file.close()


# write calculated potentials to file
def save_potentials_to_file(active_site_points: list, results_folder: str, file_name: str):

    potentials_folder = os.path.join(results_folder, 'potentials')
    if not os.path.exists(potentials_folder):
        os.makedirs(potentials_folder)

    with open(os.path.join(results_folder, f'{file_name[:-4]}.csv'), 'w') as file:
        file.write('x,y,z,protein_coulomb,protein_lennard_jones,ligand_coulomb,ligand_lennard_jones\n')
        for point in active_site_points:
            file.write(
                f'{point["coordinates"][0]},{point["coordinates"][1]},{point["coordinates"][2]},'
                f'{point["protein_coulomb"]},{point["protein_lennard_jones"]},'
                f'{point["ligand_coulomb"]},{point["ligand_lennard_jones"]}'
                '\n')
        file.close()