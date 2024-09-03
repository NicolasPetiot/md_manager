from .md_files import *
from .parameters import ATOMIC_MASSES, ATOM_NAME_SELECTION_CHI

import pandas as pd
import numpy as np

import sys
from urllib.request import urlopen
from os import path

__all__ = ["load", "pdb2df", "fetch_PDB", "atomic_masses", "com", "shift_df", "rotate_df", "angles", "dihedral_angles", "chain_theta_angles", "chain_gamma_angles", "residue_chi_angles"]

######################################################################################################################################
# Load structure

def load(filename:str) -> pd.DataFrame:
    """
    Returns a DataFrame associated to the atoms found in the first frame of the input file.
    """
    file_format = {
        ".xyz" : XYZ,
        ".pdb" : PDB,
        ".gro" : GRO
    }
    _, ext = path.splitext(filename)
    if not ext in file_format:
        raise ValueError(f"Unknown extention '{ext}'. File must be in format {list(file_format.keys())}")

    # No matter the file format, getting the next df must always be as follows:
    traj = file_format[ext](filename) # Creates traj
    traj = iter(traj)                 # Open file
    return next(traj)                 # Read coordinates

def pdb2df(filename, atom_only = False) -> pd.DataFrame:
    """
    Retruns a DataFrame associated to the atoms found in the first model of the input PDB file.
    """
    pdb = PDB(filename)
    pdb = iter(pdb)
    df = next(pdb)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df

def fetch_PDB(pdb_code:str, atom_only = False) -> pd.DataFrame:
    """
    Returns a pandas.DataFrame associated to the structure loaded from the 'rcsb.org' website.
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    df = PDB.next_frame(lines)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df

######################################################################################################################################
# DataFrame manipulation

def atomic_masses(df:pd.DataFrame) -> pd.Series:
    """Returns a Series containing atomic masses in Dalton unit based on the 'e' (element) columns of the input df."""
    def get_atom_mass(elem:str) -> float:
        try:
            m = ATOMIC_MASSES[elem.upper()]

        except KeyError:
            m = 1.0
        return m
    
    if not "e" in df.columns:
        df["e"] = ""

    elem = df["e"]
    mass = elem.apply(get_atom_mass)

    num_unknown_elem = (mass == 1.0).sum()
    if num_unknown_elem > 0:
        print(f"Warning: Found {num_unknown_elem}/{len(mass)} unknown elements. The input DataFrame should contains 'e' (element) column. The known elements are {list(ATOMIC_MASSES.keys())}")

    return mass

def com(df:pd.DataFrame) -> pd.DataFrame:
    """Returns the center of mass of the DataFrame"""
    xyz = ["x", "y", "z"]
    pos = df[xyz]

    if "m" in df.columns:
        mass = df["m"]
    else:
        mass = atomic_masses(df)

    return pos.apply(lambda x, m: x*m, m=mass).sum() / mass.sum()    

def shift_df(df:pd.DataFrame, translation_vector:np.ndarray|pd.Series) -> pd.DataFrame:
    """Shift the coordinates xyz of the DataFarme by a translation vector"""
    # Get position and apply shifting using vectorized methods
    xyz = ["x", "y", "z"]
    pos = df[xyz]
    df[xyz] = pos.apply(lambda s, v: s+v, v=translation_vector, axis = 1)

def rotate_df(df:pd.DataFrame, axis:np.ndarray|pd.Series, angle:float) -> pd.DataFrame:
    """Rotate the coordinates xyz of the DataFrame around an axis and of an angle"""
    #TODO
    raise NotImplemented("Must implement general rotation operation")

######################################################################################################################################
# Conformational angles

class InvalidChainException(Exception):
    """Raised when a chain is missing / contains too many atoms"""
    pass

def check_chain_validity(chain:pd.DataFrame, maxwarn = 5) -> None:
    """"""
    will_raise_exception = False
    if not "resi" in chain:
        print("Warning: Input DataFrame does not contains 'resi' column. Will assume that the associated chain is not missing residue.")

    else:
        Nwarn = 0
        for win in chain.resi.rolling(2):
            if len(win) > 1:
                resi1 = win.iloc[1]
                resi0 = win.iloc[0]
                diff  = resi1 - resi0
                
                if diff == 0:
                    print(f"Warning: Found several atoms that bellongs to the same residue at index {resi0}.")
                    will_raise_exception = True
                    Nwarn += 1
                
                if diff < 0:
                    print("Warning: Found several chains in the input DataFrame.")
                    will_raise_exception = True
                    Nwarn += 1
                
                if diff > 1:
                    print(f"Warning: Missing residue(s) between resi {resi0:4d} and {resi1:4d}")
                    will_raise_exception = True
                    Nwarn += 1

                if Nwarn == maxwarn:
                    break

    if will_raise_exception:
        raise InvalidChainException("Found warning(s) that does not allows conformational angles calculation.")

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

def chain_theta_angles(chain:pd.DataFrame, maxwarn = 10) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Theta angles of the input chain. The series index is the `resi` number 
    and the associated value is the computed theta angles formed by the surounding "CA" atoms.

    Can be used as `frame.groupby(...)[["name", "resi", "x", "y", "z"]].apply(chain_theta_angles)` for multimers.
    """
    check_chain_validity(chain, maxwarn)

    xyz = ["x", "y", "z"]
    pos = chain[xyz].to_numpy()

    theta = pd.Series(index=chain.index, dtype=float, name="Theta")
    theta.iloc[1:-1] = angles(pos)

    return theta

def chain_gamma_angles(chain:pd.DataFrame, maxwarn = 10) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Gamma angles of the input chain. The series index is the `resi` number 
    and the associated value is the computed gamma angles formed by the surounding "CA" atoms.

    Can be used as `frame.groupby(...)[["name", "resi", "x", "y", "z"]].apply(chain_gamma_angles)` for multimers.
    """
    check_chain_validity(chain, maxwarn)

    xyz = ["x", "y", "z"]
    pos = chain[xyz].to_numpy()

    gamma = pd.Series(index=chain.index, dtype=float, name="Gamma")
    gamma.iloc[1:-2] = dihedral_angles(pos)

    return gamma

def residue_chi_angles(res:pd.DataFrame) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Chi angles of the input residue. The series index in the 'ChiN' identifier 
    and the associated value is the computed chi angle formed by the surounding atoms.

    Can be used as `frame.groupby(...)[["name", "resn", "x", "y", "z"]].apply(residue_chi_angles)`.
    """
    xyz = ["x", "y", "z"]
    atm_selection = ATOM_NAME_SELECTION_CHI[res.resn.unique()[0]]
    res = res.query("name in @atm_selection")
    chis = dihedral_angles(res[xyz].to_numpy())
    chis = [chi for chi in chis] + [np.nan for _ in range(5 - len(chis))]
    
    return pd.Series(chis, index=["Chi%d"%i for i in range(1, 6)])

