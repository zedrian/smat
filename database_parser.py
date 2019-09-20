from pandas import read_csv
from classes import ResidueDesc, AtomDesc, ResiduesDatabase
import os


# construct database from the parametric files
# todo make DB file ones to speed up the work and only fill the residues with atoms
def get_residues_description() -> ResiduesDatabase:
    # get atom types description ones for all other functions if needed
    # some aa have 2 variants of protonation state
    # if 'res_f_p' (residues for print) are given it writes csv file for visualisation

    # get some data from csv table
    data = read_csv('Database/Atoms_properties.csv', header=0, delimiter=';')

    # initiate database class object
    residues_database = ResiduesDatabase()

    # construct ResiduDesc objects from the prep file
    def construct_resdesclist_from_prep(residues_database: ResiduesDatabase, prepfile: str):

        # divide prep file on described residues per lines
        elements = dict()
        index = 0
        elements[index] = []
        with open(prepfile, 'r') as file:
            for line in file.readlines()[2:]:
                if line.rstrip() != 'DONE' and line.rstrip() != 'STOP':
                    elements[index].append(line.rstrip())
                elif line.rstrip() == 'DONE':
                    index += 1
                    elements[index] = []
                    continue
                else:
                    break
            file.close()

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
            if elements[key][2].split(' ')[0] != '':
                short_name = elements[key][2].split(' ')[0]
                long_name = None
            else:
                short_name = elements[key][2].split(' ')[1]
                long_name = elements[key][0]

            res_desc = ResidueDesc(long_name=long_name, short_name=short_name, atoms=list())
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

                        # short the number of atom types with the common properties
                        atom_type = ''
                        for atom in range(0, len(data)):
                            # fill the properties for common atoms
                            if data['Atom'][atom] == atom_desc.get_type().upper():
                                atom_type = list(data['Atom']).index(atom_desc.get_type().upper())
                        if atom_desc.get_type().upper()[0] == 'C' and atom_desc.get_type().upper() != 'CZ':
                            atom_type = list(data['Atom']).index('C*')
                        elif atom_desc.get_type().upper()[0] == 'N':
                            atom_type = list(data['Atom']).index('N')
                        elif atom_desc.get_type().upper()[0] == 'H':
                            atom_type = list(data['Atom']).index('HO')
                        elif atom_desc.get_type().upper()[0] == 'O':
                            atom_type = list(data['Atom']).index('O2')
                        elif atom_desc.get_type().upper()[0] == 'P':
                            atom_type = list(data['Atom']).index('P')

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

            residues_database.add_residue(short_name=res_desc.get_short_name(), residue=res_desc)

    # construct ResiduDesc objects from the libraries
    def construct_resdesclist_from_lib(residues_database: ResiduesDatabase):
        libdir = 'Database/lib'
        for root, dirs, filenames in os.walk(libdir):
            for file in filenames:
                lines = open(os.path.join(root, file), 'r').readlines()
                res_name = lines[1].rstrip()[1:]
                res_desc = ResidueDesc(long_name=res_name.replace('"', ''), short_name=res_name.replace('"', ''), atoms=list())
                res_atoms_lines = []
                for line in lines[3:]:
                    if line[0] != '!':
                        res_atoms_lines.append(line.rstrip())
                    else:
                        break
                for line in res_atoms_lines:
                    line_elements = [x.strip('\"') for x in line.split(' ')[1:]]
                    atom_desc = AtomDesc(fullname=line_elements[0], type=line_elements[1], parent_name=res_desc.get_short_name(), charge=line_elements[7])

                    # short the number of atom types with the common properties
                    atom_type = ''
                    for atom in range(0, len(data)):
                        # fill the properties for common atoms
                        if data['Atom'][atom] == atom_desc.get_type().upper():
                            atom_type = list(data['Atom']).index(atom_desc.get_type().upper())
                    if atom_desc.get_type().upper()[0] == 'C' and atom_desc.get_type().upper() != 'CZ':
                        atom_type = list(data['Atom']).index('C*')
                    elif atom_desc.get_type().upper()[0] == 'N':
                        atom_type = list(data['Atom']).index('N')
                    elif atom_desc.get_type().upper()[0] == 'H':
                        atom_type = list(data['Atom']).index('HO')
                    elif atom_desc.get_type().upper()[0] == 'O':
                        atom_type = list(data['Atom']).index('O2')
                    elif atom_desc.get_type().upper()[0] == 'P':
                        atom_type = list(data['Atom']).index('P')

                    atom_desc.radius = float(data['Radius'][atom_type])
                    atom_desc.edep = float(data['Edep'][atom_type])
                    atom_desc.A = float(data['A'][atom_type])
                    atom_desc.B = float(data['B'][atom_type])
                    atom_desc.mass = float(data['Mass'][atom_type])

                    # fill the properties for atoms that are not normal
                    # we have undefined atom type CZ - for sp hybridisated C

                    res_desc.atoms.append(atom_desc)

                residues_database.add_residue(short_name=res_desc.get_short_name(), residue=res_desc)

    prepdir = 'Database/prep/'
    for root, dirs, filenames in os.walk(prepdir):
        for prepfile in filenames:
            construct_resdesclist_from_prep(residues_database, os.path.join(root, prepfile))
    construct_resdesclist_from_lib(residues_database)

    # copy the properties of HIS to another form of the HIS
    residues_database.add_residue(short_name='HIS', residue=residues_database.get_residue('HIE'))
    # residues_database.add_residue(short_name='HIP', residue=residues_database.get_residue('HIS'))

    return residues_database


database = get_residues_description()
