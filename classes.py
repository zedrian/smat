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


class ResiduesDatabase:
    def __init__(self, residues: dict = dict(), physical_atoms = list()):
        self.residues = residues
        self.amino_acids = self.construct_amino_acids_list()
        self.cofactors = self.construct_cofactors_list()
        self.physical_atoms = physical_atoms

    def add_residue(self, short_name: str, residue):
        self.residues[short_name] = residue

    def get_residue(self, residue_short_name: str):
        return self.residues[residue_short_name]

    def construct_amino_acids_list(self) -> list:
        amino_acids = list()

        for key in self.residues:
            if self.residues[key].get_amino_acid_letter() is not None:
                amino_acids.append(key)

        return amino_acids

    def construct_cofactors_list(self) -> list:
        cofactors = list()

        for key in self.residues:
            if self.residues[key].cofactor:
                cofactors.append(key)

        return cofactors

    def get_amino_acids(self):
        return self.amino_acids

    def get_cofactors(self):
        return self.cofactors

    def get_physical_atoms(self):
        return self.physical_atoms

    def get_physical_atom(self, pattern):  # pattern is either id or coords (list)
        if isinstance(pattern, UUID):
            for atom in self.physical_atoms:
                if atom.get_id() == pattern:
                    return atom
        elif isinstance(pattern, list):
            for atom in self.physical_atoms:
                if atom.get_coords() == pattern:
                    return atom


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

    def get_fullname(self):
        return self.fullname

    def get_type(self):
        return self.type

    def get_radius(self):
        return self.radius

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

    def get_ids(self):
        return self.ids

    def add_id(self, id: UUID):
        self.ids.append(id)


class PhysicalAtom:
    def __init__(self, bio_atom: Atom, atom_desc: AtomDesc, coords = list()):
        self.bio_atom = bio_atom
        self.atom_desc = atom_desc
        self.coords = coords
        self.id = uuid4()

    def __repr__(self):
        return f'Atom coords: {self.coords}\n' + \
                self.atom_desc.__repr__()

    def get_x(self):
        return self.coords[0]

    def get_y(self):
        return self.coords[1]

    def get_z(self):
        return self.coords[2]

    def get_coords(self):
        return self.coords

    def get_atom_desc(self):
        return self.atom_desc

    def get_id(self):
        return self.id


class ResidueDesc:
    def __init__(self, short_name: str, atoms: list, long_name: str = None):
        self.long_name = long_name
        self.short_name = short_name
        self.atoms = atoms
        self.amino_acid_letter = self.if_amino_acid(self.short_name)
        self.cofactor = self.if_cofactor(self.short_name)

    def __repr__(self):
        return f'Residue name: {self.long_name}\n' + \
            f'Residue short name: {self.short_name}\n' + \
            f'Residue atoms: {[x.get_type() for x in self.atoms]}'

    def get_long_name(self):
        return self.long_name

    def get_short_name(self):
        return self.short_name

    def get_atoms(self) -> list:
        return self.atoms

    def get_atom(self, atom_fullname: str) -> AtomDesc:
        for atom in self.atoms:
            if atom.get_fullname() == atom_fullname:
                return atom

        print(f'ERROR: attempt to get atom {atom_fullname} from residue {self.long_name}')
        return None

    def get_amino_acid_letter(self):
        return self.amino_acid_letter

    def if_cofactor(self):
        return self.cofactor

    # define if the residue is amino acid
    def if_amino_acid(self, short_name: str) -> str:
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

        if short_name != '':
            try:
                return amino_acids[short_name]
            except KeyError:
                return None
        else:
            return None

    # define if the residue is cofactor
    def if_cofactor(self, short_name: str) -> bool:
        cofactors = ['GDP', 'GTP', 'ADP', 'ATP', 'FMN', 'FAD', 'NAD', 'HEM']

        if short_name in cofactors:
            return True
        else:
            return False