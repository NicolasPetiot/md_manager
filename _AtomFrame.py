from ._parameters import *
import pandas as pd
import numpy as np

__all__ = ["MdTraj", "pdb2df", "df2pdb", "atom_position", "atom_mass", "atom_bfactors"]
class MdTraj :
    """
    Iterator object that allow to loop over the frames in a md trajectory in the PDB format.
    """
    def __init__(self, filename) -> None:
        self.file = open(filename, "r")
        self.file.close()

    def __len__(self) -> int :
        """
        Returns the number of frames in trajectory
        """
        len = 0 
        if self.file.closed :
            self.file = open(self.file.name, "r")
            for line in self.file :
                if line[:5] == "MODEL":
                    len += 1
            self.file.close()

        else :
            for line in self.file :
                if line[:5] == "MODEL":
                    len += 1

        return len

    def __enter__(self):
        """
        Called using 'with MdTraj(filename) as self'.
        """
        self.file = open(self.file.name, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Called at the end of a 'with MdTraj(filename) as self' statment.
        """
        self.file.close()

    def __iter__(self):
        """
        Return the iterator object itself.
        """
        if self.file.closed:
            self.file = open(self.file.name, "r")
        return self
    
    def __next__(self):
        """
        Return the next frame
        """
        atoms = []
        for line in self.file :
            if line[:6] == "ENDMDL":
                return pd.DataFrame(atoms, columns=["name", "resn", "chain", "resi", "segi", "alt", "x", "y", "z", "b", "e", "m", "hetatm"])
            
            # append atoms with potential new informations
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                atoms.append(scan_pdb_line(line))
            
        # end of file :
        self.file.close()
        raise StopIteration
    
    def __del__(self):
        """
        Called with del self.
        """
        self.file.close()


def pdb2df(filename:str) -> pd.DataFrame:
    """
    Returns a pandas.DataFrame that contains informations about all the atoms in a pdb file.
    """
    with open(filename, "r") as file :
        atoms = []
        for line in file :
            if line[:4] == "ATOM" or line[:6] == "HETATM" :
                atoms.append(scan_pdb_line(line))

    return pd.DataFrame(atoms, columns=["name", "resn", "chain", "resi", "segi", "alt", "x", "y", "z", "b", "e", "m", "hetatm"])

def df2pdb(df:pd.DataFrame, filename:str):
    """
    Generates a pdb file from the atoms in inputed DataFrame.
    """
    with open(filename, "w") as file :
        file.write("HEADER Generated with md_manager.")
        file.write("MODEL     1\n")
        
        id = 0
        APO = df.query("hetatm == False")
        HET = df.query("hetatm == True")

        for _, chain in APO.groupby(["chain"]):
            for _, atom in chain.iterrows():
                id += 1
                file.write(serie2pdb_line(atom, id))
            file.write("TER\n")
        
        for _, atom in HET.iterrows():
            id += 1
            file.write(serie2pdb_line(atom, id))
        file.write("ENDMDL\n")


def serie2pdb_line(s:pd.Series, id:int):
    """
    Generates a string that contains informations about an atom in pdb format.
    """

    if s.hetatm :
        hetatm = "HETATM"
    else :
        hetatm = "ATOM  "

    # strings :
    name = s["name"]
    alt = s["alt"]
    resn = s["resn"]
    chain = s["chain"]
    resi = s["resi"]
    segi = s["segi"]
    e = s["e"]
    # insertion code missing

    # floating numbers
    x = s.x
    y = s.y
    z = s.z
    # occupancy missing
    b = s.b

    return f"{hetatm}{id:5d} {name}{alt}{resn} {chain}{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}      {b:6.2f}      {segi}{e}\n"


def scan_pdb_line(line):
    """
    Returns a tuple that contains the (name, resn, chain, resi, segi, alt, x, y, z, beta, symbol, mass, hetatm) of the line in pdf file.
    If the line does not start with ATOM or HETATM, the method raise a ValueError.
    If the atom element is not included in the ATOMIC_MASSES parameter, the method raise a KeyError.
    """
    if line[:4] == "ATOM":
        hetatm = False
    elif line[:6] == "HETATM":
        hetatm = True
    else :
        raise ValueError("line does not contains atom's information.")
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
    try :
        mass = ATOMIC_MASSES[symbol]
    except KeyError :
        raise KeyError(f"Unknown atomic symbol {symbol}, please update the ATOMIC_MASSES dictionary.")

    return name, resn, chain, resi, segi, alt, x, y, z, beta, symbol, mass, hetatm


def atom_position(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the [x, y, z] coordinates of the atoms in df.
    """
    return df[["x", "y", "z"]].to_numpy(dtype=float)

def atom_mass(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the masses of the atoms in df.
    """
    return df["m"].to_numpy(dtype=float)

def atom_bfactors(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the thermal bfactors of the atoms in df.
    """
    return df["b"].to_numpy(dtype=float)

    