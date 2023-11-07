from ._parameters import *
import numpy as np
from copy import deepcopy

class Atom :

    def __init__(self, chain = " ", resn = "   ", resi = 0, name = "    ", x = 0.0, y = 0.0, z = 0.0, hetatm = False, beta = 0.0, segi = "    ", alt = " ", elem = " ") -> None:
        """
        The Atom object contains the informations related to atoms in files such as pdb. It aims at providing easy acess and manipulation of such files.
        """
        self.chain = chain
        self.resn  = resn
        self.resi  = resi
        self.name  = name

        self.x = x
        self.y = y
        self.z = z
        self.xyz = np.array([x, y, z])

        self.hetatm = hetatm

        self.e = elem
        self.m = ATOMIC_MASSES[self.e]
        self.b = beta

        self.segi = segi
        self.alt  = alt # alternate location indicator

    def __repr__(self) -> str:
        """
        Executed with print(self).
        """
        return f"{self.chain}-{self.resn[0].upper()}{self.resn[1:].lower()}{self.resi:<3} {self.name} [{self.x:8.3f}; {self.y:8.3f}; {self.z:8.3f}]"
    
def line2atom(line:str):
    """
    Read a line of a pdb file and create the associated Atom object.
    """
    #atom_num = int(line[6:11])
    name     = line[12:16]
    alt      = line[16] 
    resn     = line[17:20]
    chain    = line[21]
    resi     = int(line[22:26])
    x        = float(line[30:38])
    y        = float(line[38:46])
    z        = float(line[46:54])
    try :
        beta = float(line[60:66])
    except ValueError :
        beta = 0.0
    segi     = line[72:76]
    symbol   = line[76:78]

    atom = Atom(chain, resn, resi, name, x, y, z, beta = beta, segi=segi, alt = alt, elem=symbol)

    return atom


class AtomSelection :
    ### Special Methods ###
    def __init__(self, atoms:[Atom] = []) -> None:
        """
        An AtomSelection object is basically a list of Atom objects with the ability to select them based on their properties.
        """
        self._atoms = []

        # generator :
        for atom in atoms :
            self.append(atom)

    def __repr__(self)-> str :
        """
        Called with print(self).
        """
        if len(self) == 0 :
            return "[]"
        out = "[\n"
        for atom in self :
            out += atom.__repr__()
            out += "\n"
        out += "]"
        return out
    
    def __len__(self) -> int:
        """
        Called with 'len(self)'.
        """
        return len(self._atoms)
    
    def __getitem__(self, index):
        """
        Called with 'self[index]'.
        """
        return self._atoms[index]
    
    def __iadd__(self, other):
        """
        Called with 'self += other'.
        """
        for atom in other :
            self.append(atom)
        return self
    
    def __add__(self, other):
        """
        Called with 'self + other'.
        """
        sele = deepcopy(self)
        sele += other

        return sele

    def append(self, atm:Atom) -> None:
        """
        Add an Atom object at the end of the AtomSelection object.
        """
        self._atoms.append(atm)

    ### selector Methods ###
    def chain(self, selector = "", **kwargs):
        """
        Returns an AtomSelection object that contains atoms that have been selected by chains.
        Args :
            selector (str)       : Expected chain for the output AtomSelection.
            not (str)            : Chain that will not be in the output AtomSelection.
            or (list[str] / str) : List of chains that are expected in the output AtomSelection.

        """
        if "Not" in kwargs :
            return AtomSelection([atm for atm in self if atm.chain != kwargs["Not"]])
        elif "Or" in kwargs :
            return AtomSelection([atm for atm in self if atm.chain == selector or atm.chain in kwargs["Or"]])
        return AtomSelection([atm for atm in self if atm.chain == selector])
    
    def resi(self, selector = "", **kwargs):
        """
        Returns an AtomSelection object that contains atoms that have been selected by residue index.
        Args :
            selector (str)       : Expected resi for the output AtomSelection.
            not (str)            : Residue index that will not be in the output AtomSelection.
            or (list[str] / str) : List of residue indexes that are expected in the output AtomSelection.

        """
        if "Not" in kwargs :
            return AtomSelection([atm for atm in self if atm.resi != kwargs["Not"]])
        elif "Or" in kwargs :
            return AtomSelection([atm for atm in self if atm.resi == selector or atm.resi in kwargs["Or"]])
        return AtomSelection([atm for atm in self if atm.resi == selector])
    
    def resn(self, selector = "", **kwargs):
        """
        Returns an AtomSelection object that contains atoms that have been selected by residue name.
        Args :
            selector (str)       : Expected resn for the output AtomSelection.
            not (str)            : Residue name that will not be in the output AtomSelection.
            or (list[str] / str) : List of residue names that are expected in the output AtomSelection.

        """
        if "Not" in kwargs :
            return AtomSelection([atm for atm in self if atm.resn != kwargs["Not"]])
        elif "Or" in kwargs :
            return AtomSelection([atm for atm in self if atm.resn == selector or atm.resn in kwargs["Or"]])
        return AtomSelection([atm for atm in self if atm.resn == selector])
    
    def name(self, selector = "", **kwargs):
        """
        Returns an AtomSelection object that contains atoms that have been selected by atom names.
        Args :
            selector (str)       : Expected name for the output AtomSelection.
            not (str)            : Atom name that will not be in the output AtomSelection.
            or (list[str] / str) : List of atom names that are expected in the output AtomSelection.

        """
        if "Not" in kwargs :
            return AtomSelection([atm for atm in self if atm.name != kwargs["Not"]])
        elif "Or" in kwargs :
            return AtomSelection([atm for atm in self if atm.name == selector or atm.name in kwargs["Or"]])
        return AtomSelection([atm for atm in self if atm.name == selector])
    
    def segi(self, selector = "", **kwargs):
        """
        Returns an AtomSelection object that contains atoms that have been selected by segment identifiers.
        Args :
            selector (str)       : Expected segi for the output AtomSelection.
            not (str)            : Segment identifier that will not be in the output AtomSelection.
            or (list[str] / str) : List of segment identifiers that are expected in the output AtomSelection.

        """
        if "Not" in kwargs :
            return AtomSelection([atm for atm in self if atm.segi != kwargs["Not"]])
        elif "Or" in kwargs :
            return AtomSelection([atm for atm in self if atm.segi == selector or atm.segi in kwargs["Or"]])
        return AtomSelection([atm for atm in self if atm.segi == selector])
    
    def hetatm(self, selector:bool):
        """
        Returns an AtomSelection object that contains atoms that have been selected by nature (atoms or hetero atoms).
        """
        return AtomSelection([atm for atm in self if atm.hetatm == selector])
    
    ### Analyse Methods ###
    def center_of_mass(self):
        """
        Returns the center of mass as well of the total mass of the atoms within the AtomSelection object.
        """
        center_of_mass = np.zeros(3)
        total_mass = 0.0

        for atom in self :
            center_of_mass += atom.xyz
            total_mass += atom._mass

        return center_of_mass/total_mass, total_mass
    
    def get_xyz(self):
        """
        Returns a numpy.ndarray that contains the atom's position in the AtomSelection object.
        """
        return np.array([atm.xyz for atm in self])
    
    def get_m(self):
        """
        Returns a numpy.ndarray that contains the atom's mass in the AtomSelection object.
        """
        return np.array([atm.m for atm in self])
    
    def get_b(self):
        """
        Returns a numpy.ndarray that contains the atom's thermal factors in the AtomSelection object.
        """
        return np.array([atm.b for atm in self])
    
    def get_chain(self):
        """
        Returns a numpy.ndarray that contains the chains name in the AtomSelection object.
        """
        return np.unique([atm.chain for atm in self])
    
    def get_resi(self):
        """
        Returns a numpy.ndarray that contains the residue indexes in the AtomSelection object.
        """
        return np.unique([atm.resi for atm in self])
    
    def get_resn(self):
        """
        Returns a numpy.ndarray that contains the residue names in the AtomSelection object.
        """
        return np.unique([atm.resn for atm in self])
    
    def get_name(self):
        """
        Returns a numpy.ndarray that contains the atom names in the AtomSelection object.
        """
        return np.unique([atm.name for atm in self])
    
    def get_segi(self):
        """
        Returns a numpy.ndarray that contains the segment identifiers in the AtomSelection object.
        """
        return np.unique([atm.segi for atm in self])

### load and write pdb files ###
def pdb2sele(filename:str):
    """
    Returns an AtomSelection object that contains all the atoms of a pdb file
    """
    sele = AtomSelection()
    with open(filename, "r") as file :
        for line in file :
            if line[:4] == "ATOM":
                atom = line2atom(line)
                sele.append(atom)
            
            elif line[:6] == "HETATM":
                atom = line2atom(line)
                atom._hetatm = True
                sele.append(atom)

    return sele

def sele2pdb(sele:AtomSelection, filename = "sele.pdb"):
    """
    Geneate a pdb file associated to the AtomSelection object.
    """
    with open(filename, "w") as file :
        file.write("REMARK generated with sele2pdb\n")
        file.write("MODEL     1\n")
        
        id = 0
        
        atoms = sele.hetatm(False)
        for s in atoms.get_chains():
            chain = atoms.chain(s)
            for atom in chain :
                id += 1
                file.write(f"ATOM  {id:5d} {atom._name:^4s}{atom._alt:1s}{atom._resn:3s} {atom._chain:1s}{atom._resi:4d}    {atom._x:8.3f}{atom._y:8.3f}{atom._z:8.3f}      {atom._beta:6.2f}                    {atom._elem:>2s}\n")
            file.write("TER\n")
            
        hetatms = sele.hetatm(True)
        for atom in hetatms :
            file.write(f"ATOM  {id:5d} {atom._name:^4s}{atom._alt:1s}{atom._resn:3s} {atom._chain:1s}{atom._resi:4d}    {atom._x:8.3f}{atom._y:8.3f}{atom._z:8.3f}      {atom._beta:6.2f}                    {atom._elem:>2s}\n")
        file.write("ENDMDL")

def scan_traj(filename: str, func, **kwargs):
    """
    Create an AtomSelection object for each step of the trajectory file. Apply a specified function to the selection for all frame in the trajectory
    """
    with open(filename, "r") as file :
        for line in file :
            if line[:5] == "MODEL":
                model = AtomSelection()

            elif line[:6] == "ENDMDL":
                func(model, **kwargs)

            elif line[:4] == "ATOM":
                atm = line2atom(line)
                model.append(atm)

            elif line[:6] == "HETATM":
                atm = line2atom(line)
                atm.hetatm = True
                model.append(atm)