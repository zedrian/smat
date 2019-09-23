from uuid import uuid4, UUID
from Bio.PDB import Atom


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
    def __init__(self, fullname: str = '', type: str = '', radius=float(), charge=float(), edep=float(), A=float(), B=float(), mass=float(), parent_name: str = '', ids=list()):
        self.fullname = fullname.replace(' ', '')
        self.type = type
        self.radius = radius
        self.charge = charge
        self.edep = edep
        self.A = A
        self.B = B
        self.mass = mass
        self.parent_name = parent_name
        self.ids = ids

    def __repr__(self):
        return f'Atom full name: {self.fullname}\n' + \
            f'Atom type: {self.type}\n' + \
            f'Atom radius: {self.radius}\n' + \
            f'Atom charge: {self.charge}\n' + \
            f'Atom edep: {self.edep}\n' + \
            f'Atom Van der Waals 12-part (A): {self.A}\n' + \
            f'Atom Van der Waals 6-part (B): {self.B}\n' + \
            f'Atom mass: {self.mass}\n' + \
            f'Atom parent: {self.parent_name}'

    def get_fullname(self) -> str:
        return self.fullname

    def get_type(self) -> str:
        return self.type

    def get_radius(self) -> float:
        return self.radius

    def get_charge(self) -> float:
        return self.charge

    def get_edep(self) -> float:
        return self.edep

    def get_A(self) -> float:
        return self.A

    def get_B(self) -> float:
        return self.B

    def get_mass(self) -> float:
        return self.mass

    def get_parent_name(self) -> str:
        return self.parent_name

    def get_ids(self) -> list: # ids of Physical Atoms
        return self.ids

    def add_id(self, id: UUID):
        self.ids.append(id)


class PhysicalAtom:
    def __init__(self, bio_atom: Atom, atom_desc: AtomDesc, coords: list):
        self.bio_atom = bio_atom
        self.atom_desc = atom_desc
        self.coords = coords
        self.id = uuid4()

    def __repr__(self):
        return f'Atom coords: {self.coords}\n' + \
                self.atom_desc.__repr__()

    def get_x(self) -> float:
        return self.coords[0]

    def get_y(self) -> float:
        return self.coords[1]

    def get_z(self) -> float:
        return self.coords[2]

    def get_coords(self) -> list:
        return self.coords

    def get_atom_desc(self) -> AtomDesc:
        return self.atom_desc

    def get_id(self) -> UUID:
        return self.id


class ResidueDesc:
    def __init__(self, short_name: str, atoms: list, long_name: str = None, ids=list(), terminus=None):
        self.long_name = long_name
        self.short_name = short_name
        self.atoms = atoms
        self.amino_acid_letter = self.if_amino_acid()
        self.cofactor = self.if_cofactor()
        self.ids = ids
        self.terminus = terminus

    def __repr__(self):
        return f'Residue name: {self.long_name}\n' + \
            f'Residue short name: {self.short_name}\n' + \
            f'Residue terminus: {self.terminus}\n' + \
            f'Residue atoms: {[x.get_type() for x in self.atoms]}'

    def get_long_name(self) -> str:
        return self.long_name

    def get_short_name(self) -> str:
        return self.short_name

    # list of atoms descriptions
    def get_atoms(self) -> list:
        return self.atoms

    # get atom description object
    def get_atom(self, atom_fullname: str) -> AtomDesc:
        for atom in self.atoms:
            if atom.get_fullname() == atom_fullname:
                return atom

        print(f'ERROR: attempt to get atom {atom_fullname} from residue {self.long_name}')
        return None

    def get_amino_acid_letter(self) -> str:
        return self.amino_acid_letter

    # define if the residue is amino acid
    def if_amino_acid(self) -> str:
        amino_acids ={'ALA': 'A',
                      'ARG': 'R',
                      'ASN': 'N',
                      'ASP': 'D',
                      'CYS': 'C',
                      'GLN': 'Q',
                      'GLU': 'E',
                      'GLY': 'G',
                      'HIS': 'H',
                      'HIE': 'H',
                      'HIP': 'H',
                      'HFD': 'H',
                      'ILE': 'I',
                      'LEU': 'L',
                      'LYS': 'K',
                      'MET': 'M',
                      'PHE': 'F',
                      'PRO': 'P',
                      'SER': 'S',
                      'THR': 'T',
                      'TRP': 'W',
                      'TYR': 'Y',
                      'VAL': 'V',
                      'CYM': 'C'}

        if self.short_name != '':
            try:
                return amino_acids[self.short_name]
            except KeyError:
                return None
        else:
            return None

    # define if the residue is cofactor
    def if_cofactor(self) -> bool:
        cofactors = ['GDP', 'GTP', 'ADP', 'ATP', 'FMN', 'FAD', 'NAD', 'HEM']

        if self.short_name in cofactors:
            return True
        else:
            return False

    def get_ids(self) -> list:
        return self.ids

    def add_id(self, id):
        self.ids.append(id)

    def if_terminus(self):
        return self.terminus

    def set_short_name(self, short_name: str):
        self.short_name = short_name
        return self


class PhysicalResidue:
    def __init__(self, index: int, residue_desc: ResidueDesc, atoms=list()):
        self.index = index
        self.residue_desc = residue_desc
        self.atoms = atoms
        self.id = uuid4()

    def get_atoms(self) -> list:  # get list of physical atoms
        return self.atoms

    def get_index(self) -> int:
        return self.index

    def get_id(self) -> UUID:
        return self.id

    def get_residue_desc(self) -> ResidueDesc:
        return self.residue_desc

    def get_atom(self, id: UUID) -> PhysicalAtom:
        for atom in self.get_atoms():
            if atom.get_id() == id:
                return atom


class ResiduesDatabase:
    def __init__(self, residues=list()):
        self.residues = residues
        self.amino_acids = self.construct_amino_acids_list()
        self.cofactors = self.construct_cofactors_list()

    def __repr__(self):
        return f'Residues: {self.residues}\n'

    def add_residue(self, residue: ResidueDesc):
        self.residues.append(residue)

    def get_residue(self, residue_short_name: str, terminus=None) -> ResidueDesc:
        for residue in self.residues:
            if residue.get_short_name() == residue_short_name and residue.if_terminus() == terminus:
                return residue

    def construct_amino_acids_list(self) -> list:
        amino_acids = list()

        for residue in self.residues:
            if residue.get_amino_acid_letter() is not None:
                amino_acids.append(residue.get_short_name())

        return amino_acids

    def construct_cofactors_list(self) -> list:
        cofactors = list()

        for residue in self.residues:
            if residue.cofactor:
                cofactors.append(residue.get_short_name())

        return cofactors

    def get_amino_acids(self) -> list:
        return self.amino_acids

    def get_cofactors(self) -> list:
        return self.cofactors