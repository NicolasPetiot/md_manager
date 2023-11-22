from ._parameters import *
import pandas as pd
import numpy as np

import sys
from urllib.request import urlopen

__all__ = ["PDBfile", "atom_position", "atom_mass", "atom_bfactors", "fetch_protein_data_bank", "angles", "dihedral_angles", "theta_angles", "gamma_angles", "chi_angles"]

class PDBfile :
    """
    Iterator object that allow to loop over the frames in a md trajectory in the PDB format.

    PDBfile(filename).read2df() returns a frame that contains all the atoms of the file.
    """
    def __init__(self, filename) -> None:
        self.file = open(filename, "r")
        self.file.close()

    def __len__(self) -> int :
        """
        Returns the number of frames in trajectory.
        """
        len = 0 
        close = False
        if self.file.closed :
            self.file = open(self.file.name, "r")
            close = True

        for line in self.file :
            if line[:5] == "MODEL":
                len += 1

        if close :
            self.file.close()

        return len

    ### Allow to use 'with statments' ###
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

    
    #####################################

    ### Allow to loop over frames ###
    def __iter__(self):
        """
        Is called at the beggining of a 'for model in pdb' statment.

        Caution : pdb.file must be oppened.
        """
        if self.file.closed:
            self.file = open(self.file.name, "r")
        return self
    
    def __next__(self):
        """
        Is called using 'for model in pdb'.
        Returns a pandas.DataFrame that corresponds to the next MODEL in PDB file.

        Caution : pdb.file must be oppened.

        Caution : the text readed must contains 'ENDMDL' at the end of a model. If not, please use 'read2df' instead.
        """
        atoms = []
        for line in self.file :
            if line[:6] == "ENDMDL":
                return _build_df_from_atom_list(atoms)
            
            # append atoms with potential new informations
            if line[:6] in {"ATOM  ", "HETATM"}:
                atoms.append(_scan_pdb_line(line))
            
        raise StopIteration
    
    def open(self):
        """
        Opens the pdb.file wrapper.
        """
        self.file = open(self.file.name, "r")

    def close(self):
        """
        Closes the pdb.file wrapper.
        """
        self.file.close()
    #################################

    ### Manipulate pdb lines ###
    def read2df(self):
        """
        Returns a pandas.DataFrame that contains all the atoms fount in the PDB file.
        """
        atoms = []
        with open(self.file.name, "r") as file :
            for line in file:
                if line[:6] in {"ATOM  ", "HETATM"}:
                    atoms.append(_scan_pdb_line(line))
        return _build_df_from_atom_list(atoms) 
    ############################
    

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

### Load from Protein Data Bank ###

def fetch_protein_data_bank(pdb_code:str)->pd.DataFrame:
    """
    Returns a pandas.DataFrame associated to the structure loaded from the 'rcsb.org' website.
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    atoms = []
    for line in lines :
        if line[:6] in {"ATOM  ", "HETATM"}:
            atoms.append(_scan_pdb_line(line))

    return _build_df_from_atom_list(atoms)


### Conformational angles computations ###

def angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Returns the angles associated to the ensemble of positions in 'atom_position'.
    The length of the output is len(atom_position) - 2 because angles are not defined for the first and last atom of a chain.
    """
    B1 = atom_position[1:-1, :] - atom_position[:-2, :] 
    B2 = atom_position[2:, :] - atom_position[1:-1, :] 

    B1 = B1 / np.sqrt(np.sum(B1**2, axis = 1, keepdims=True))
    B2 = B2 / np.sqrt(np.sum(B2**2, axis = 1, keepdims=True))

    return np.rad2deg(np.arccos(np.sum(-B1 * B2, axis=1)))

def dihedral_angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Returns the dihedral angles associated to the ensemble of positions in 'atom_position'.
    The length of the output is len(atom_position) - 3 because dihedral angles are not defined for the first, last and second last atom of a chain.
    """
    B1 = atom_position[1:-2, :] - atom_position[:-3, :] 
    B2 = atom_position[2:-1, :] - atom_position[1:-2, :] 
    B3 = atom_position[3:, :] - atom_position[2:-1, :]

    B1 = B1 / np.sqrt(np.sum(B1**2, axis = 1, keepdims=True))
    B2 = B2 / np.sqrt(np.sum(B2**2, axis = 1, keepdims=True))
    B3 = B3 / np.sqrt(np.sum(B3**2, axis = 1, keepdims=True))

    N1 = np.cross(B1, B2)
    N2 = np.cross(B2, B3)

    cos = np.sum(N1 * N2, axis = 1)
    sin = np.sum(np.cross(N1, N2)*B2, axis = 1)

    return np.rad2deg(np.arctan2(sin, cos))

def theta_angles(df:pd.DataFrame):
    """
    Returns a pandas.Series that contains theta angles associated to the chains in inputed df.
    """
    theta = []
    for _, chain in df.query("name == ' CA '").groupby(["chain"]):
        theta.append(np.nan)
        xyz = chain[["x", "y", "z"]].to_numpy(dtype=float)
        theta += list(angles(xyz))
        theta.append(np.nan)
    return pd.Series(theta, index=df.index, name="theta")

def gamma_angles(df:pd.DataFrame):
    """
    Returns a pandas.Series that contains the gamma angles associated to the chains in inputed df.
    """
    gamma = []
    for _, chain in df.query("name == ' CA '").groupby(["chain"]):
        gamma.append(np.nan)
        xyz = atom_position(chain)
        gamma += list(dihedral_angles(xyz))
        gamma.extend([np.nan, np.nan])
    return pd.Series(gamma, index=df.index, name="gamma")

def chi_angles(df:pd.DataFrame):
    """
    Returns a dictionary that contains the chi angles associated to the residues in inputed df.
    """
    Chis = {"%d"%i : [] for i in range(1, 6)}
    for _, residue in df.groupby(["chain", "resi"]):
        # select only atoms for the computation of chi angles
        resn = residue["resn"].unique()[0]
        residue = residue.query("name in @ATOM_NAME_SELECTION[@resn]") 

        # chi angle computation :
        xyz = atom_position(residue)
        chis = dihedral_angles(xyz)
        for i, chi in enumerate(chis):
            i += 1
            Chis["%d"%i].append(chi)
        for i in range(len(chis) + 1, 6):
            Chis["%d"%i].append(np.nan)

    return Chis

### Useful functions that will not be in the __all__ list ###

def _build_df_from_atom_list(atoms:list):
    """
    Returns a pandas.DataFrame that contains the atoms in the input list.
    """
    columns=["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q", "m"]
    return pd.DataFrame(atoms, columns=columns)

def _scan_pdb_line(line:str):
    """
    Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) 
    informations about an atom line in the PDB file.
    """
    if line.startswith("ATOM"):
        record_name = "ATOM  "
    elif line.startswith("HETATM"):
        record_name = "HETATM"
    else :
        raise ValueError("Input 'line' is not associated to an atom in a pdb file.")
    
    #atom_id   = line[ 6:11]
    name      = line[12:16]
    alt       = line[16:17]
    resn      = line[17:20]
    chain     = line[21:22]
    resi      = int(line[22:26])
    insertion = line[26:27]
    x         = float(line[30:38])
    y         = float(line[38:46])
    z         = float(line[46:54])
    try :
        occupancy = float(line[54:60])
    except ValueError :
        occupancy = 1.0
    try :
        b         = float(line[60:66])
    except ValueError :
        b         = 0.0
    segi      = line[72:76]
    elem      = line[76:78]
    charge    = line[78:79]
    try :
        mass      = ATOMIC_MASSES[elem]
    except KeyError:
        raise KeyError(f"Unknown element symbol '{elem}'. Please update the ATOMIC_MASSES dictionary in '_parameters.py'.")
    return record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass

def _generate_pdb_line(atom:pd.Series, id:int) -> str:
    """
    Generates a pdb line from an input series that contains atom's information.
    """
    # extract strings :
    record_name = atom["record_name"]
    name        = atom["name"]
    alt         = atom["alt"]
    resn        = atom["resn"]
    chain       = atom["chain"]
    insertion   = atom["insertion"]
    segi        = atom["segi"]
    elem        = atom["e"]
    charge      = atom["q"]

    # extract float/int
    id    = f"{id:5d}"
    resi  = f"{atom.resi:4d}"
    x     = f"{atom.x:8.3f}"
    y     = f"{atom.y:8.3f}"
    z     = f"{atom.z:8.3f}"
    occup = f"{atom.occupancy:6.2f}"
    b     = f"{atom.b:6.2f}"

    line = "%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"%(
        record_name, id, name, alt, resn, chain, resi, insertion, x, y, z, occup, b, segi, elem, charge
    )
    return line

def df2pdb(cls, df:pd.DataFrame, filename: str, title = "", remark = "") -> None:
    """
    Generates a pdb file from an inputed DataFrame.
    """
    with open(filename, "w") as file :
        if title != "":
            file.write("TITLE     " + title + "\n")
        if remark != "":
            for line in remark.splitlines():
                file.write("REMARK    " + line + "\n")
            
        file.write("MODEL        1\n")
        id = 0
        APO = df[df.record_name == "ATOM  "]
        HET = df[df.record_name == "HETATM"]
        for _, chain in APO.groupby(["chain"]):
            for _, atom in chain.iterrows():
                id += 1
                line = cls._generate_pdb_line(atom, id)
                file.write(line)
            file.write("TER\n")
        for atom in HET.iterrows():
            id += 1
            line = _generate_pdb_line(atom, id)
            file.write(line)
        file.write("TER\n")
        file.write("ENDMDL\n")