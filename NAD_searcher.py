from Bio.PDB import *
from Bio.PDB.Chain import Chain
from pandas import read_csv, DataFrame
from os import path
import time


# return a chain from main data file
def get_chain(file_name: str, chain: str) -> Chain:

    parser = PDBParser(QUIET=True)

    # get structure
    structure = parser.get_structure('pdb', file_name)

    # get chain
    chain = structure[0][chain]
    return chain


# return a list of atoms that belong to pattern sequence
def get_pattern_atoms(residues: list, pattern_sequence: str, amino_acid_code: DataFrame) -> list:
    # construct one-letter sequence
    sequence = ''
    for residue in residues:
        for index in range(len(amino_acid_code['three'])):
            if residue.get_resname().lower() == amino_acid_code['three'].tolist()[index]:
                sequence += amino_acid_code['one'][index]

    # find index of pattern coords
    pattern_start = residues[sequence.find(pattern_sequence)].get_full_id()[-1][1]
    pattern_end = pattern_start + len(pattern_sequence)

    # get all atoms of the pattern
    atoms = []

    # get all atoms from residues of pattern
    for residue in residues:
        if residue.get_full_id()[-1][1] in range(pattern_start, pattern_end):
            res_atoms = [a for a in residue.get_atoms()]
            for a in res_atoms:
                atoms.append(a)

    return atoms


# return a list of ligand (I hope only NAD(P)) atoms
def get_ligand_atoms(residues: list, amino_acid_code: DataFrame) -> list:
    for residue in residues:
        if not residue.get_resname().lower() in amino_acid_code['three'].tolist() and residue.get_resname().lower() != 'coa':
            atoms = [a for a in residue.get_atoms() if not a.get_fullname().startswith('H')]
            if 48 >= len(atoms) >= 41:
                return atoms
    return list()


# return atoms of contact
def get_contacts(pattern_atoms: list, ligand_atoms: list) -> list:
    ns = NeighborSearch(pattern_atoms)

    contacts = list()
    for atom in ligand_atoms:
        close_atoms = ns.search(atom.coord, 4)
        if len(close_atoms) > 0:
            for a in close_atoms:
                if not a in contacts:
                    contacts.append(a)

    return contacts

start_time = time.time()

# get single-letter amino acid code
amino_acid_code = read_csv('aas.csv', ',', header=0)

folder_path = r'D:\work\filtered-pdb'
pdb_hits = read_csv('pdb-hits.filtered.csv', header=0)
total_contacts = list()

for row_index, row in pdb_hits.iterrows():
    file_path = path.join(folder_path, f'{row["PDB_code"].lower()}.pdb')
    chain = row['Chain']
    pattern_sequence = row['Pattern_sequence']

    if path.isfile(file_path):
        # get residues without water
        chain = get_chain(file_path, chain)
        residues = [r for r in chain.get_residues() if r.get_resname() != 'HOH']

        # run functions
        ligand_atoms = get_ligand_atoms(residues, amino_acid_code)
        all_other_atoms = [a for a in chain.get_atoms() if a not in ligand_atoms]
        contacts = get_contacts(all_other_atoms, ligand_atoms)
        total_contacts.append(', '.join([f'{atom.name}[{atom.get_parent().get_id()[1]}]:{atom.coord}' for atom in contacts]))
        print(f'\r{row_index+1} / {pdb_hits.shape[0]}', end='')
    else:
        print(f'\r{row["PDB_code"]} not found')
        total_contacts.append('FILE_NOT_FOUND')
print()

# store contacts in dataframe
pdb_hits['Contacts'] = total_contacts
# remove rows without contacts
pdb_hits.drop(pdb_hits[pdb_hits['Contacts'] == ''].index, inplace=True)

pdb_hits.to_csv('contacts.4.csv')
print(f'done in {time.time() - start_time} seconds')
