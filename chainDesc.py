from Bio.PDB import Chain
from classes import ResidueDesc, PhysicalAtom
from database_parser import database

class ChainDesc:
    def __init__(self, chain: Chain, ligand: ResidueDesc = None):
        self.chain = chain
        self.ligand = ligand
        self.residues = self.get_transformed_residues()
        self.has_ligands = self.contains_ligands()
        self.short_seq = self.get_shortened_sequence()

    def get_transformed_residues(self):
        transformed_residues = list()

        for residue in self.chain.get_residues():
            atoms_desc = list()

            for atom in residue.get_atoms():
                atom_desc = database.get_residue(atom.get_parent().get_resname()).get_atom(atom.get_id())  # get atom description for physical atom from database
                physical_atom = PhysicalAtom(coords=atom.get_coord(), atom_desc=atom_desc, bio_atom=atom)  # generate object
                database.physical_atoms.append(physical_atom)  # add physical atoms to database
                atom_desc.add_id(physical_atom.get_id())
                atoms_desc.append(atom_desc)

            res_desc = ResidueDesc(short_name=residue.get_resname(), atoms=atoms_desc)
            transformed_residues.append(res_desc)

        return transformed_residues

    def get_residues(self):
        return self.residues

    def contains_ligands(self) -> bool:
        for residue in self.residues:
            if not residue.get_short_name() in database.get_amino_acids() + database.get_cofactors():
                atoms = [a for a in residue.get_atoms() if not a.get_name().startswith('H')]
                if len(atoms) >= 6:
                    return True
        return False

    def get_shortened_sequence(self) -> str:
        line = ''
        for residue in self.chain.get_residues():
            if database.get_residue(residue.get_resname()).get_amino_acid_letter() is None:
                continue
            line += database.get_residue(residue.get_resname()).get_amino_acid_letter()
        return line

    def get_chain_length(self) -> int:
        return len(self.short_seq)

    def if_has_ligands(self) -> bool:
        return self.has_ligands

    def get_short_seq(self) -> str:
        return self.short_seq

    def get_bio_chain(self):
        return self.chain