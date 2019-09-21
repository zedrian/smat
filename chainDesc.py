from uuid import UUID
from Bio.PDB import Chain
from classes import PhysicalAtom, PhysicalResidue
from database_parser import database


class ChainDesc:
    def __init__(self, chain: Chain, ligand: PhysicalResidue = None):
        self.chain = chain
        self.ligand = ligand
        self.residues = self.get_transformed_residues()
        self.has_ligands = self.contains_ligands()
        self.shortened_seq = self.get_shortened_sequence()

    def get_transformed_residues(self) -> list:   # list of Physical Residues
        transformed_residues = list()

        for residue in self.chain.get_residues():  # Biopython Residue
            index = residue.get_segid()
            terminus = None
            if index == 1:
                terminus = 'N'
            elif 'OXT' in [a.get_id() for a in residue.get_atoms()]:
                terminus = 'C'

            residue_desc = database.get_residue(residue.get_resname(), terminus=terminus)
            physical_residue = PhysicalResidue(index=index, residue_desc=residue_desc)
            physical_atoms = list()
            for atom in residue.get_atoms():  # Biopython Atom!
                atom_desc = residue_desc.get_atom(atom_fullname=atom.get_id())  # get atom description for physical atom from database
                physical_atom = PhysicalAtom(coords=atom.get_coord(), atom_desc=atom_desc, bio_atom=atom)  # generate object
                atom_desc.add_id(physical_atom.get_id())
                physical_atoms.append(physical_atom)

            physical_residue.atoms = physical_atoms
            residue_desc.add_id(physical_residue.get_id())  # add physical residue id to residue description list of ids

            transformed_residues.append(physical_residue)

        return transformed_residues

    def get_residues(self) -> list:  # list of Physical Residues
        return self.residues

    def get_residue_by_id(self, id: UUID) -> PhysicalResidue:
        for residue in self.residues:
            if residue.get_id() == id:
                return residue

    def contains_ligands(self) -> bool:
        for residue in self.residues:
            if not residue.get_residue_desc().get_short_name() in database.get_amino_acids() + database.get_cofactors():
                atoms = [a for a in residue.get_atoms() if not a.get_atom_desc().get_fullname().startswith('H')]
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
        return len(self.shortened_seq)

    def if_has_ligands(self) -> bool:
        return self.has_ligands

    def get_shortened_seq(self) -> str:
        return self.shortened_seq

    def get_bio_chain(self) -> Chain:
        return self.chain