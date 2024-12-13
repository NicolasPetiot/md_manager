from .md_files import *
from .parameters import ATOMIC_MASSES

import pandas as pd
import numpy as np

import sys
from urllib.request import urlopen
from os import path

__all__ = ["load", "save", "pdb2df", "fetch_PDB", "atomic_masses", "com", "shift_df", "rotate_df"]

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

def save(filename:str, data:pd.DataFrame|list[pd.DataFrame]):
    """
    Creates trajectory file from input data.
    """
    file_format = {
        ".xyz" : XYZ,
        ".pdb" : PDB,
        ".gro" : GRO
    }
    _, ext = path.splitext(filename)
    if not ext in file_format:
        raise ValueError(f"Unknown extention '{ext}'. File must be in format {list(file_format.keys())}")
    
    new = file_format[ext](filename, "w")
    if type(data) == pd.DataFrame:
        data = [data]

    for i, df in enumerate(data):
        new.write_frame(df, model_id=i+1)


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

def df2pdb(filename:str, data:pd.DataFrame|list[pd.DataFrame]):
    """Generates a pdb structure/trajctory according to the type of the input data"""
    new = PDB(filename, "w")

    if type(data) == pd.DataFrame:
        data = [data]

    for i, df in enumerate(data):
        new.write_frame(df, model_id=i+1)

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
                    print("Maximum warning reached...")

    if will_raise_exception:
        raise InvalidChainException("Found warning(s) that does not allows conformational angles calculation.")